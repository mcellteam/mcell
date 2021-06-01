/******************************************************************************
 *
 * Copyright (C) 2021 by
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

#ifndef SRC4_REGION_EXPR_H_
#define SRC4_REGION_EXPR_H_

#include "defines.h"

namespace MCell {

class World;

enum class RegionExprOperator {
  INVALID,
  UNION,
  INTERSECT,
  DIFFERENCE,
  LEAF
};


class RegionExprNode {
public:
  RegionExprNode()
    : op(RegionExprOperator::INVALID),
      region_id(REGION_INDEX_INVALID),
      left(nullptr), right(nullptr) {
  }

  ~RegionExprNode() {
    // children are contained in ReleaseEvent::all_region_expr_nodes,
    // and are deleted when ReleaseEvent is destroyed
  }

  RegionExprOperator op;

  region_id_t region_id;

  RegionExprNode* left;
  RegionExprNode* right;

  void dump(const World* world) const; // does not print any newlines
  std::string to_string(const World* world, const bool for_datamodel = false) const;
};


class RegionExpr {
public:
  RegionExpr() :
    root(nullptr) {
  }
  ~RegionExpr();

  // constructor and container for all region expr nodes
  RegionExprNode* create_new_region_expr_node_leaf(const region_id_t region_id);
  RegionExprNode* create_new_region_expr_node_op(const RegionExprOperator op, RegionExprNode* left, RegionExprNode* right);

  // used when release_shape is ReleaseShape::Region
  RegionExprNode* root;

private:
  std::vector<RegionExprNode*> all_region_expr_nodes;
};

} /* namespace MCell */

#endif /* SRC4_REGION_EXPR_H_ */
