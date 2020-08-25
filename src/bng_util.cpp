/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#include "bng_util.h"

#include <vector>
#include <cassert>
#include <map>

using namespace std;

struct Node {
  bool is_mol; // not component
  string name;
  string state; // for component
  string compartment; // for molecule
  vector<size_t> connections;
};


// used in set to represent edges
// copied from rxn_rule.xpp
class UnorderedPair {
public:
  UnorderedPair(const int a, const int b) : first(std::min(a,b)), second(std::max(a,b)) {
  }
  bool operator == (const UnorderedPair& b) const {
    return first == b.first && second == b.second;
  }
  bool operator < (const UnorderedPair& b) const {
    if (first < b.first) {
      return true;
    }
    else if (first == b.first) {
      return second < b.second;
    }
    else {
      return false;
    }
  }

  int first;
  int second;
};


static Node convert_node_str(const string& s) {
  assert(s.size() > 0);
  Node res;

  if (s[0] == 'c') {
    res.is_mol = false;
  }
  else if (s[0] == 'm') {
    res.is_mol = true;
  }
  else {
    assert(false);
  }

  if (!res.is_mol) {
    // component

    size_t tilde = s.find('~');
    assert(tilde != string::npos);
    res.name = s.substr(2, tilde - 2);

    size_t excl = s.find('!');
    assert(excl != string::npos); // component must be connected
    res.state = s.substr(tilde + 1, excl - (tilde + 1));
    if (res.state == "NO_STATE") {
      res.state = "";
    }
  }
  else {
    // molecule
    size_t at = s.find('@');
    res.name = s.substr(2, at - 2);

    size_t excl = s.find('!');
    if (excl == string::npos) {
      excl = s.size();
    }

    res.compartment = s.substr(at + 1, excl - (at + 1));
    if (res.compartment == "NOCOMPARTMENT") {
      res.compartment = "";
    }
  }

  // connections
  size_t excl = s.find('!');
  if (excl != string::npos) {
    size_t pos = excl + 1;
    do {
      excl = s.find('!', pos);
      if (excl == string::npos) {
        excl = s.size();
      }

      string num = s.substr(pos, excl - pos);
      res.connections.push_back(stoi(num));

      pos = excl + 1;
    } while (excl != s.size());
  }

  return res;
}


static void graph_pattern_to_nodes(const string pat, vector<Node>& nodes) {

  size_t pos = 0;
  size_t comma_or_end;
  do {

    comma_or_end = pat.find(',', pos);
    if (comma_or_end == string::npos) {
      comma_or_end = pat.size();
    }

    if (comma_or_end == string::npos) {
      comma_or_end = pat.size();
    }

    string node_str = pat.substr(pos, comma_or_end - pos);
    Node n = convert_node_str(node_str);
    nodes.push_back(n);

    pos = comma_or_end;
    pos++;
  } while (comma_or_end != pat.size() && pos != pat.size());
}


static std::string nodes_to_bngl(const vector<Node>& nodes) {
  string res;

  map<UnorderedPair, int> bonds;

  for (size_t mol_index = 0; mol_index < nodes.size(); mol_index++) {
    const Node& mol = nodes[mol_index];
    if (mol.is_mol) {
      res += mol.name;
      if (!mol.connections.empty()) {
        res += "(";

        for (size_t i = 0; i < mol.connections.size(); i++) {
          size_t component_index = mol.connections[i];
          const Node& comp = nodes[component_index];

          res += comp.name;
          if (comp.state != "") {
            res += "~" + comp.state;
          }

          for (size_t conn: comp.connections) {
            if (conn == mol_index) {
              // connection to myself
              continue;
            }
            UnorderedPair bond(component_index, conn);
            int bond_value;
            auto it = bonds.find(bond);
            if (it != bonds.end()) {
              bond_value = it->second;
            }
            else {
              int next_bond_index = bonds.size() + 1; // we start from 1
              bonds[bond] = next_bond_index;
              bond_value = next_bond_index;
            }
            res += "!" + to_string(bond_value);
          }

          if (i + 1!= mol.connections.size()) {
            res += ",";
          }
        }
        res += ")";
      }

      if (mol.compartment != "") {
        res += "@" + mol.compartment;
      }

      // expecting that all the molecules are at the end
      if (mol_index + 1 != nodes.size()) {
        res += ".";
      }
    }
  }

  return res;
}


std::string graph_pattern_to_bngl(const char* graph_pattern) {
  assert(graph_pattern != nullptr);

  // each c: or m: is a node on graph, ! specifies to what it is connected
  // c:a~R!2!1,c:b~T!3!0,m:A@NOCOMPARTMENT!0,m:B@NOCOMPARTMENT!1

  vector<Node> nodes;

  // process the nauty string
  graph_pattern_to_nodes(graph_pattern, nodes);

  // make BNGL representation
  return nodes_to_bngl(nodes);
}

