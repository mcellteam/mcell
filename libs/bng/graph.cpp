/*
 * graph.cpp
 *
 *  Created on: Jun 22, 2020
 *      Author: ahusar
 */

#include "graph.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

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
      const Graph& graph1_, const Graph& graph2_, MultipleMappingsVector& mappings_)
    : graph1(graph1_), graph2(graph2_), mappings(mappings_) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) {

    mappings.push_back(MappingVector());
    BGL_FORALL_VERTICES_T(v, graph1, Graph) {
      mappings.back().push_back(
          std::make_pair(
              get(boost::vertex_index_t(), graph1, v) ,
              get(boost::vertex_index_t(), graph2, get(f, v))));
    }

    return true;
  }
  std::vector<MappingVector> get_setvmap() { return mappings; }

private:
  const Graph& graph1;
  const Graph& graph2;
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
void get_subgraph_isomorphism_mappings(Graph& pattern, Graph& cplx, MultipleMappingsVector& res) {

  res.clear();

  // setting result to store the resulting mappings
  CallBackToCollectMapping callback(pattern, cplx, res);

  PropertyMapMoleculeTypeMatching vertex_comp =
      make_property_map_molecule_type_matching(get(vertex_name, pattern), get(vertex_name, cplx));

  boost::vf2_subgraph_iso(
      pattern, cplx, std::ref(callback),
      boost::vertex_order_by_mult(pattern),
      boost::vertices_equivalent(vertex_comp)
  );

#ifdef DEBUG_CPLX_MATCHING
  int i = 0;
  for (std::vector <std::pair<graph_type::vertex_descriptor, graph_type::vertex_descriptor>>& vec: callback.get_setvmap()) {
    cout << i << ": ";
    for (std::pair<graph_type::vertex_descriptor, graph_type::vertex_descriptor> item: vec) {
      cout << "(" << item.first << "," << item.second << ")" << " ";
    }

    i++;
    cout << "\n";
  }
#endif
}

} // namespace BNG
