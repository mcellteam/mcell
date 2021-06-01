/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#include "region_expr.h"
#include "datamodel_defines.h"
#include "world.h"

using namespace std;

namespace MCell {

string RegionExprNode::to_string(const World* world, const bool for_datamodel) const {
  stringstream out;
  assert(op != RegionExprOperator::INVALID);

  if (op == RegionExprOperator::LEAF) {
    const string& region_name = world->get_region(region_id).name;
    if (for_datamodel) {
      return DMUtils::get_object_w_region_name(region_name);
    }
    else {
      return region_name;
    }
  }

  assert(left != nullptr);
  out << "(";
  out << left->to_string(world, for_datamodel);

  switch(op) {
    case RegionExprOperator::UNION:
      out << " + ";
      break;
    case RegionExprOperator::INTERSECT:
      out << " * ";
      break;
    case RegionExprOperator::DIFFERENCE:
      out << " - ";
      break;
    default:
      assert(false);
  }
  out << right->to_string(world, for_datamodel);
  out << ")";
  return out.str();
}


void RegionExprNode::dump(const World* world) const {
  cout << to_string(world);
}


RegionExprNode* RegionExpr::create_new_region_expr_node_leaf(const region_id_t region_id) {
  RegionExprNode* res = new RegionExprNode;
  res->op = RegionExprOperator::LEAF;
  res->region_id = region_id;
  all_region_expr_nodes.push_back(res);
  return res;
}


RegionExprNode* RegionExpr::create_new_region_expr_node_op(
    const RegionExprOperator op, RegionExprNode* left, RegionExprNode* right) {

  RegionExprNode* res = new RegionExprNode;
  res->op = op;
  res->left = left;
  res->right = right;
  all_region_expr_nodes.push_back(res);
  return res;
}


RegionExpr::~RegionExpr() {
  for (RegionExprNode* expr_node: all_region_expr_nodes) {
    delete expr_node;
  }
}

} /* namespace MCell */
