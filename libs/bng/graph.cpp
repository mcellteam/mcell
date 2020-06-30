/*
 * graph.cpp
 *
 *  Created on: Jun 22, 2020
 *      Author: ahusar
 */

#include "graph.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "debug_config.h"

using namespace std;
using namespace boost;

namespace BNG {

std::ostream & operator<<(std::ostream &out, const Node& n) {
  if (n.is_mol) {
    std::cout << "m:" << n.mol->mol_type_id;
  }
  else {
    cout << "c:" << n.component->component_type_id;
    if (n.component->state_id == STATE_ID_DONT_CARE) {
      cout << "~" << "DONT_CARE";
    }
    else {
      cout << "~" << n.component->state_id;
    }
    if (n.component->bond_value == BOND_VALUE_ANY) {
      cout << "!+";
    }
    else if (n.component->bond_value == BOND_VALUE_NO_BOND) {
      cout << "!NO_BOND";
    }
    else {
      cout << "!" << n.component->bond_value;
    }
  }
  return out;
}



class CallBackToCollectMapping {

public:
  typedef Graph::vertex_descriptor vertex_descriptor;

  // constructor, result is stored into mappings_
  CallBackToCollectMapping(
      const Graph& graph1_, const Graph& graph2_, const bool only_first_match_, MultipleMappingsVector& mappings_)
    : graph1(graph1_), graph2(graph2_), only_first_match(only_first_match_), mappings(mappings_) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) {

    // TODO: handle maximal number of matches, some counter

    mappings.push_back(MappingVector());
    BGL_FORALL_VERTICES_T(v, graph1, Graph) {
      mappings.back()[get(boost::vertex_index_t(), graph1, v)] =
              get(boost::vertex_index_t(), graph2, get(f, v));
    }

    if (only_first_match) {
      return false; // terminates search
    }
    else {
      return true;
    }
  }
  std::vector<MappingVector> get_setvmap() { return mappings; }

private:
  const Graph& graph1;
  const Graph& graph2;
  bool only_first_match;
  MultipleMappingsVector& mappings; // result is stored here
};

// Binary function object that returns true if the values for item1
// in property_map1 and item2 in property_map2 are equivalent.
// TODO: add references where possible
struct PropertyMapMoleculeTypeMatching {

  PropertyMapMoleculeTypeMatching(const VertexNameMap& property_map1_, const VertexNameMap& property_map2_) :
    property_map1(property_map1_),
    property_map2(property_map2_) { }

  template <typename ItemFirst, typename ItemSecond>
  bool operator()(const ItemFirst item1, const ItemSecond item2) {
    const Node& n1 = get(property_map1, item1);
    const Node& n2 = get(property_map2, item2);
#ifdef DEBUG_CPLX_MATCHING
    cout << "Comparing " << item1 << ": " << n1 << " and " << item2 << ": " << n2 << "\n";
#endif
    bool res = n1.compare(n2);
#ifdef DEBUG_CPLX_MATCHING
    cout << "  " << (res ? "true" : "false") << "\n";
#endif
    return res;
  }

private:
  const VertexNameMap& property_map1;
  const VertexNameMap& property_map2;
};


// Returns a property_map_equivalent object that compares the values
// of property_map1 and property_map2.
PropertyMapMoleculeTypeMatching make_property_map_molecule_type_matching(
    const VertexNameMap& property_map1,
    const VertexNameMap& property_map2
) {
  return PropertyMapMoleculeTypeMatching(property_map1, property_map2);
}


// TODO LATER: can we make the arguments constant?
// using mutable graph now
void get_subgraph_isomorphism_mappings(Graph& pattern, Graph& cplx, const bool only_first_match, MultipleMappingsVector& res) {

  res.clear();

  // setting result to store the resulting mappings
  CallBackToCollectMapping callback(pattern, cplx,  only_first_match, res);

  PropertyMapMoleculeTypeMatching vertex_comp =
      make_property_map_molecule_type_matching(get(vertex_name, pattern), get(vertex_name, cplx));

  boost::vf2_subgraph_iso(
      pattern, cplx, std::ref(callback),
      boost::vertex_order_by_mult(pattern),
      boost::vertices_equivalent(vertex_comp)
  );

#ifdef DEBUG_CPLX_MATCHING
  int i = 0;
  for (MappingVector& vec: callback.get_setvmap()) {
    cout << i << ": ";
    for (std::pair<Graph::vertex_descriptor, Graph::vertex_descriptor> item: vec) {
      cout << "(" << item.first << "," << item.second << ")" << " ";
    }

    i++;
    cout << "\n";
  }
#endif
}

void dump_graph(const Graph& g_const) {

  // not manipulating with the graph in any way, but boost needs
  // non-const variant
  Graph& g = const_cast<Graph&>(g_const);
  VertexNameMap index = boost::get(boost::vertex_name, g);

  // for each molecule instance in pattern_graph
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
  pair<vertex_iter, vertex_iter> it;
  for (it = boost::vertices(g); it.first != it.second; ++it.first) {

    Graph::vertex_descriptor desc = *it.first;
    const Node& mol = index[desc];
    if (mol.is_mol) {
      cout << (int)desc << ": " << mol << "\n";

      // also print connected components
      graph_traits<Graph>::out_edge_iterator ei, edge_end;
      for (boost::tie(ei,edge_end) = boost::out_edges(desc, g); ei != edge_end; ++ei) {
        Graph::edge_descriptor e_mol_comp = *ei;
        Graph::vertex_descriptor comp_desc = boost::target(e_mol_comp, g);

        const Node& comp = index[comp_desc];
        assert(!comp.is_mol && "Only a component may be connected to a molecule.");
        cout << " -> " << (int)comp_desc << ": " << comp << ", connections: ";

        // and to which components they are connected
        graph_traits<Graph>::out_edge_iterator ei_comp, edge_end2;
        for (boost::tie(ei_comp,edge_end2) = boost::out_edges(comp_desc, g); ei_comp != edge_end2; ++ei_comp) {
          Graph::edge_descriptor e_comp_comp = *ei_comp;
          Graph::vertex_descriptor connected_comp_desc = boost::target(e_comp_comp, g);

          const Node& comp = index[connected_comp_desc];
          if (!comp.is_mol) {
            cout << connected_comp_desc << ", ";
          }
        }
        cout << "\n";
      }
    }
  }
}

} // namespace BNG
