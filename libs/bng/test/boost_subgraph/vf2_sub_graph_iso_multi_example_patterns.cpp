
// based on boost example vf2_sub_graph_iso_multi_example.cpp
//
// g++ -std=c++14 vf2_sub_graph_iso_multi_example_patterns.cpp -O0 -g3 -o example.elf

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
using namespace boost;
using namespace std;


namespace Graph {

// molecule instance/pattern is a node in the complex graph, edge is a bond
// note: char is replaced with a pointer to MoleculeInstance in the actual implementation
typedef property< vertex_name_t, char, property< vertex_index_t, int > > mt_vertex_property;

// using a vecS graphs => the index maps are implicit.
// we are not using any edge property
typedef adjacency_list< vecS, vecS, bidirectionalS, mt_vertex_property> graph_type;

typedef property_map< graph_type, vertex_name_t >::type vertex_name_map_t;

// Binary function object that returns true if the values for item1
// in property_map1 and item2 in property_map2 are equivalent.
struct property_map_molecule_type_matching {

  property_map_molecule_type_matching(const vertex_name_map_t property_map1, const vertex_name_map_t property_map2) :
    m_property_map1(property_map1),
    m_property_map2(property_map2) { }

  template <typename ItemFirst, typename ItemSecond>
  bool operator()(const ItemFirst item1, const ItemSecond item2) {
    char n1 = get(m_property_map1, item1);
    char n2 = get(m_property_map2, item2);
    cout << "Comparing " << item1 << " " << n1 << " and " << item2 << " " << n2 << "\n";

    return n1 == n2;
  }

private:
  const vertex_name_map_t m_property_map1;
  const vertex_name_map_t m_property_map2;
};

// Returns a property_map_equivalent object that compares the values
// of property_map1 and property_map2.
property_map_molecule_type_matching make_property_map_molecule_type_matching(
    const vertex_name_map_t property_map1,
    const vertex_name_map_t property_map2
) {
  return property_map_molecule_type_matching(property_map1, property_map2);
}


} // namespace Graph


int main()
{
  // Build graph_small
  Graph::graph_type graph_small;


  add_vertex(Graph::mt_vertex_property('b'), graph_small); //0
  add_vertex(Graph::mt_vertex_property('c'), graph_small); //1
  add_vertex(Graph::mt_vertex_property('d'), graph_small); //2

  add_edge(0, 1, graph_small);
  add_edge(1, 2, graph_small);

  // Build graph_large
  Graph::graph_type graph_large;

  add_vertex(Graph::mt_vertex_property('a'), graph_large); //0
  add_vertex(Graph::mt_vertex_property('b'), graph_large); //1
  add_vertex(Graph::mt_vertex_property('c'), graph_large); //2
  add_vertex(Graph::mt_vertex_property('c'), graph_large); //3
  add_vertex(Graph::mt_vertex_property('d'), graph_large); //3

  add_edge(0, 1, graph_large);
  add_edge(1, 2, graph_large);
  add_edge(1, 3, graph_large);
  add_edge(2, 4, graph_large);
  add_edge(3, 4, graph_large);

  // Create callback
  vf2_print_callback< Graph::graph_type, Graph::graph_type > callback(graph_small, graph_large);

  //typedef Graph::property_map_molecule_type_matching< Graph::vertex_name_map_t, Graph::vertex_name_map_t > vertex_comp_t;
  typedef Graph::property_map_molecule_type_matching vertex_comp_t;

  vertex_comp_t vertex_comp = Graph::make_property_map_molecule_type_matching(get(vertex_name, graph_small), get(vertex_name, graph_large));

  // Print out all subgraph isomorphism mappings between graph1 and graph_small.
  // Function vertex_order_by_mult is used to compute the order of
  // vertices of graph1. This is the order in which the vertices are examined
  // during the matching process.
  // might need to use mono

  vf2_subgraph_iso(
      graph_small, graph_large, callback, vertex_order_by_mult(graph_small),
      vertices_equivalent(vertex_comp)
  );

  return 0;
}
