
// based on boost example vf2_sub_graph_iso_multi_example.cpp
//
// g++ -std=c++14 vf2_sub_graph_iso_multi_rxn_match.cpp -O0 -g3 -o example.elf

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
using namespace boost;
using namespace std;

typedef uint state_id_t;

const uint STATE_ID_DONT_CARE = UINT32_MAX - 2;

typedef uint bond_value_t;

const uint BOND_VALUE_ANY = UINT32_MAX - 2; // for '+' in patterns such as a!+
const uint BOND_VALUE_NO_BOND = UINT32_MAX - 3;

typedef uint component_type_id_t;
typedef uint mol_type_id_t;

namespace Graph {

struct Node {
  // default ctor that sets ids to invalid?

  static Node CreateMol(
      const mol_type_id_t mol_type_id_
  ) {
    Node res;
    res.is_mol = true;
    res.mol_type_id = mol_type_id_;
    return res;
  }

  static Node CreateComponent(
      const component_type_id_t component_type_id_,
      const state_id_t state_id_,
      const bond_value_t bond_value_
  ) {
    Node res;
    res.is_mol = false;
    res.component_type_id = component_type_id_;
    res.state_id = state_id_;
    res.bond_value = bond_value_;
    return res;
  }

  bool compare(const Node& n2) const {
    const Node& n1 = *this;

    if (n1.is_mol != n2.is_mol) {
      return false;
    }
    else if (n1.is_mol) {
      assert(n2.is_mol);
      // molecule
      return n1.mol_type_id == n2.mol_type_id;
    }
    else {
      assert(!n1.is_mol && !n2.is_mol);
      // component
      if (n1.component_type_id != n2.component_type_id) {
        // must be the same
        return false;
      }

      // state
      if (n1.state_id != STATE_ID_DONT_CARE && n2.state_id != STATE_ID_DONT_CARE) {
        // must be the same or don't care for one of the compared nodes
        if (n1.state_id != n2.state_id) {
          return false;
        }
      }

      // bond
      // comparing !+
      if (n1.bond_value == BOND_VALUE_ANY) {
        // it is ok when the second node has !+ as well
        return n2.bond_value != BOND_VALUE_NO_BOND;
      }
      if (n2.bond_value == BOND_VALUE_ANY) {
        return n1.bond_value != BOND_VALUE_NO_BOND;
      }

      // we do not care about actual bond values because to what is the component connected is
      // represented by the graph
      return true;
    }
  }


  bool is_mol;
  // molecule
  mol_type_id_t mol_type_id; // the only attribute for a molecule

  // component
  component_type_id_t component_type_id;
  state_id_t state_id;
  bond_value_t bond_value;
};

/*static inline*/ std::ostream & operator<<(std::ostream &out, const Node& n) {
  if (n.is_mol) {
    cout << "m:" << n.mol_type_id;
  }
  else {
    cout << "c:" << n.component_type_id;
    if (n.state_id == STATE_ID_DONT_CARE) {
      cout << "~" << "DONT_CARE";
    }
    else {
      cout << "~" << n.state_id;
    }
    if (n.bond_value == BOND_VALUE_ANY) {
      cout << "!+";
    }
    else if (n.bond_value == BOND_VALUE_NO_BOND) {
      cout << "!NO_BOND";
    }
    else {
      cout << "!" << n.bond_value;
    }
  }
  return out;
}


// molecule instance/pattern is a node in the complex graph, edge is a bond
typedef property< vertex_name_t, Node, property< vertex_index_t, int > > mt_vertex_property;

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
    const Node& n1 = get(m_property_map1, item1);
    const Node& n2 = get(m_property_map2, item2);
    cout << "Comparing " << item1 << ": " << n1 << " and " << item2 << ": " << n2 << "\n";

    bool res = n1.compare(n2);
    cout << "  " << (res ? "true" : "false") << "\n";
    return res;
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


using namespace Graph;

typedef std::vector <std::pair<graph_type::vertex_descriptor, graph_type::vertex_descriptor>> MappingVector;

class my_call_back {

public:
  typedef graph_type::vertex_descriptor vertex_descriptor;


  // constructor
  my_call_back(const graph_type& graph1, const graph_type& graph2) : graph1_(graph1), graph2_(graph2) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) {

    set_of_vertex_iso_map.push_back(MappingVector());
    BGL_FORALL_VERTICES_T(v, graph1_, graph_type) {
      set_of_vertex_iso_map.back().push_back(
          std::make_pair(
              get(boost::vertex_index_t(), graph1_, v) ,
              get(boost::vertex_index_t(), graph2_, get(f, v))));
    }

    return true;
  }
  std::vector<MappingVector>& get_setvmap() { return set_of_vertex_iso_map; }

private:
  const graph_type& graph1_;
  const graph_type& graph2_;
  std::vector<MappingVector> set_of_vertex_iso_map;
};

int main()
{
  // Build graph_small
  Graph::graph_type graph_complex; // molecule

  // A(c0~R,c1~T,c0~S)
  // using different numbers just for clarity
  add_vertex(Graph::mt_vertex_property(Node::CreateMol(0/*A*/)), graph_complex); //0
  add_vertex(Graph::mt_vertex_property(Node::CreateComponent(1/*c0*/, 2/*R*/, BOND_VALUE_NO_BOND)), graph_complex); //1
  add_vertex(Graph::mt_vertex_property(Node::CreateComponent(3/*c1*/, 4/*T*/, BOND_VALUE_NO_BOND)), graph_complex); //2
  add_vertex(Graph::mt_vertex_property(Node::CreateComponent(1/*c0*/, 5/*S*/, BOND_VALUE_NO_BOND)), graph_complex); //3

  add_edge(0, 1, graph_complex);
  add_edge(0, 2, graph_complex);
  add_edge(0, 3, graph_complex);

  // Build graph_large
  Graph::graph_type graph_pattern; // pattern

  add_vertex(Graph::mt_vertex_property(Node::CreateMol(0/*A*/)), graph_pattern); //0
  add_vertex(Graph::mt_vertex_property(Node::CreateComponent(1/*c0*/, STATE_ID_DONT_CARE, BOND_VALUE_NO_BOND)), graph_pattern); //1

  add_edge(0, 1, graph_pattern);


  my_call_back callback(graph_pattern, graph_complex);

  // Create callback
  //vf2_print_callback< Graph::graph_type, Graph::graph_type > callback(graph_pattern, graph_complex);

  //typedef Graph::property_map_molecule_type_matching< Graph::vertex_name_map_t, Graph::vertex_name_map_t > vertex_comp_t;
  typedef Graph::property_map_molecule_type_matching vertex_comp_t;

  vertex_comp_t vertex_comp =
      Graph::make_property_map_molecule_type_matching(get(vertex_name, graph_pattern), get(vertex_name, graph_complex));

  // Print out all subgraph isomorphism mappings between graph1 and graph_small.
  // Function vertex_order_by_mult is used to compute the order of
  // vertices of graph1. This is the order in which the vertices are examined
  // during the matching process.
  // might need to use mono

  vf2_subgraph_iso(
      graph_pattern, graph_complex, std::ref(callback),
      vertex_order_by_mult(graph_pattern),
      vertices_equivalent(vertex_comp)
  );

  int i = 0;
  for (std::vector <std::pair<graph_type::vertex_descriptor, graph_type::vertex_descriptor>>& vec: callback.get_setvmap()) {
    cout << i << ": ";
    for (std::pair<graph_type::vertex_descriptor, graph_type::vertex_descriptor> item: vec) {
      cout << "(" << item.first << "," << item.second << ")" << " ";
    }

    i++;
    cout << "\n";
  }

  return 0;
}
