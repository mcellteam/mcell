/*
 * graph.cpp
 *
 *  Created on: Jun 22, 2020
 *      Author: ahusar
 */

#include "graph.h"

#include <sstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "bng/bng_data.h"

#include "debug_config.h"



using namespace std;
using namespace boost;

namespace BNG {

string Node::to_str(const BNGData* bng_data) const {
  stringstream out;

  if (is_mol) {
    out << "m:";
    if (bng_data != nullptr) {
      out << bng_data->get_elem_mol_type(mol->elem_mol_type_id).name;
    }
    else {
      out << mol->elem_mol_type_id;
    }
  }
  else {
    out << "c:";
    if (bng_data != nullptr) {
      out << bng_data->get_component_type(component->component_type_id).name;
    }
    else {
      out << component->component_type_id;
    }
    if (component->state_id == STATE_ID_DONT_CARE) {
      out << "~" << "DONT_CARE";
    }
    else {
      out << "~";
      if (bng_data != nullptr) {
        out << bng_data->get_state_name(component->state_id);
      }
      else {
        out << component->state_id;
      }

    }
    if (component->bond_value == BOND_VALUE_BOUND) {
      out << "!+";
    }
    else if (component->bond_value == BOND_VALUE_ANY) {
      out << "!?";
    }
    else if (component->bond_value == BOND_VALUE_UNBOUND) {
      out << "!NO_BOND";
    }
    else {
      out << "!" << component->bond_value;
    }
  }
  return out.str();
}


ostream & operator<<(ostream &out, const Node& n) {
  out << n.to_str();
  return out;
}



class CallBackToCollectMapping {

public:
  typedef Graph::vertex_descriptor vertex_descriptor;

  // constructor, result is stored into mappings_
  CallBackToCollectMapping(
      const Graph& graph1_, const Graph& graph2_, const bool only_first_match_, VertexMappingVector& mappings_)
    : graph1(graph1_), graph2(graph2_), only_first_match(only_first_match_), mappings(mappings_) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) {

    // TODO: handle maximal number of matches, some counter

    mappings.push_back(VertexMapping());
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

private:
  const Graph& graph1;
  const Graph& graph2;
  bool only_first_match;
  VertexMappingVector& mappings; // result is stored here
};

// Binary function object that returns true if the values for item1
// in property_map1 and item2 in property_map2 are equivalent.
// TODO: add references where possible
struct PropertyMapMoleculeTypeMatching {

  PropertyMapMoleculeTypeMatching(VertexNameMap property_map1_, VertexNameMap property_map2_) :
    property_map1(property_map1_),
    property_map2(property_map2_) { }

  template <typename ItemFirst, typename ItemSecond>
  bool operator()(const ItemFirst item1, const ItemSecond item2) {
    const Node& n1 = get(property_map1, item1);
    const Node& n2 = get(property_map2, item2);
#ifdef DEBUG_CPLX_MATCHING_EXTRA_COMPARE
    cout << "Comparing " << item1 << ": " << n1 << " and " << item2 << ": " << n2 << "\n";
#endif
    bool res = n1.compare(n2);
#ifdef DEBUG_CPLX_MATCHING_EXTRA_COMPARE
    cout << "  " << (res ? "true" : "false") << "\n";
#endif
    return res;
  }

private:
  VertexNameMap property_map1;
  VertexNameMap property_map2;
};


// Returns a property_map_equivalent object that compares the values
// of property_map1 and property_map2.
PropertyMapMoleculeTypeMatching make_property_map_molecule_type_matching(
    VertexNameMap property_map1,
    VertexNameMap property_map2
) {
  return PropertyMapMoleculeTypeMatching(property_map1, property_map2);
}


// TODO LATER: can we make the arguments constant?
// using mutable graph now
void get_subgraph_isomorphism_mappings(
    Graph& pattern,
    Graph& cplx,
    const bool only_first_match,
    VertexMappingVector& res) {

  res.clear();

#ifdef DEBUG_CPLX_MATCHING
  cout << "\nPattern:\n";
  dump_graph(pattern);

  cout << "Cplx:\n";
  dump_graph(cplx);
#endif

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
  for (VertexMapping& vec: res) {
    cout << i << ": ";
    for (pair<Graph::vertex_descriptor, Graph::vertex_descriptor> item: vec) {
      cout << "(" << item.first << "," << item.second << ")" << " ";
    }

    i++;
    cout << "\n";
  }
#endif
}

// bng_data might be nullptr
void dump_graph(const Graph& g_const, const BNGData* bng_data, const std::string ind) {

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
      cout << ind << (int)desc << ": " << mol.to_str(bng_data) << "\n";

      // also print connected components
      graph_traits<Graph>::out_edge_iterator ei, edge_end;
      for (boost::tie(ei,edge_end) = boost::out_edges(desc, g); ei != edge_end; ++ei) {
        Graph::edge_descriptor e_mol_comp = *ei;
        Graph::vertex_descriptor comp_desc = boost::target(e_mol_comp, g);

        const Node& comp = index[comp_desc];
        assert(!comp.is_mol && "Only a component may be connected to a molecule.");
        cout << ind << " -> " << (int)comp_desc << ": " << comp.to_str(bng_data) << ", connections: ";

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

        if (comp.modified_ordering_index != INDEX_INVALID) {
          cout << "mod_ord: " << comp.modified_ordering_index;
        }
        cout << "\n";
      }
    }
  }
}


void dump_graph_mapping(const VertexMapping& mapping) {
  for (auto it: mapping) {
    cout << "(" << (int)it.first << ", " << it.second << "), ";
  }
}


} // namespace BNG
