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
#include "api/python_export_constants.h"

using namespace std;

namespace MCell {

string RegionExprNode::to_string(const World* world, const bool for_datamodel) const {
  stringstream out;
  assert(op != RegionExprOperator::INVALID);

  if (op == RegionExprOperator::LEAF_SURFACE_REGION) {
    if (world != nullptr) {
      const string& region_name = world->get_region(region_id).name;
      if (for_datamodel) {
        return DMUtils::get_object_w_region_name(region_name);
      }
      else {
        return region_name;
      }
    }
    else {
      return "region id: " + std::to_string(region_id);
    }
  }
  else if (op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    if (world != nullptr) {
      const string& go_name = world->get_geometry_object(geometry_object_id).name;
      if (for_datamodel) {
        return DMUtils::remove_obj_name_prefix(go_name) + API::REGION_ALL_SUFFIX;
      }
      else {
        return go_name;
      }
    }
    else {
      return "geometry object id: " + std::to_string(geometry_object_id);
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


RegionExprNode* RegionExpr::create_new_expr_node_leaf_surface_region(const region_id_t id) {
  assert(id != REGION_ID_INVALID);
  RegionExprNode* res = new RegionExprNode;
  res->op = RegionExprOperator::LEAF_SURFACE_REGION;
  res->region_id = id;
  all_region_expr_nodes.push_back(res);
  return res;
}


RegionExprNode* RegionExpr::create_new_expr_node_leaf_geometry_object(const geometry_object_id_t id) {
  assert(id != GEOMETRY_OBJECT_ID_INVALID);
  RegionExprNode* res = new RegionExprNode;
  res->op = RegionExprOperator::LEAF_GEOMETRY_OBJECT;
  res->geometry_object_id = id;
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


static RegionExprNode* clone_region_expr_recursively(
    RegionExpr* re,
    const RegionExprNode* src_node) {
  assert(src_node != nullptr);

  if (src_node->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    return re->create_new_expr_node_leaf_geometry_object(src_node->geometry_object_id);
  }
  else if (src_node->op == RegionExprOperator::LEAF_SURFACE_REGION) {
    return re->create_new_expr_node_leaf_surface_region(src_node->region_id);
  }
  else if (src_node->has_binary_op()) {
    return re->create_new_region_expr_node_op(
        src_node->op,
        clone_region_expr_recursively(re, src_node->left),
        clone_region_expr_recursively(re, src_node->right));
  }
  else {
    assert(false);
    return nullptr;
  }
}


RegionExpr& RegionExpr::operator=(const RegionExpr& other) {
  if (other.root == nullptr) {
    root = nullptr;
  }
  else {
    root = clone_region_expr_recursively(this, other.root);
  }
  return *this;
}


RegionExpr::~RegionExpr() {
  for (RegionExprNode* expr_node: all_region_expr_nodes) {
    delete expr_node;
  }
}

} /* namespace MCell */
