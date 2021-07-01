/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_REGION_EXPR_H_
#define SRC4_REGION_EXPR_H_

#include "defines.h"

namespace MCell {

class World;
class Region;

enum class RegionExprOperator {
  INVALID,
  UNION,
  INTERSECT,
  DIFFERENCE,
  LEAF_SURFACE_REGION,
  LEAF_GEOMETRY_OBJECT
};


class RegionExprNode {
public:
  RegionExprNode()
    : op(RegionExprOperator::INVALID),
      region_id(REGION_INDEX_INVALID),
      geometry_object_id(GEOMETRY_OBJECT_ID_INVALID),
      left(nullptr), right(nullptr) {
  }

  ~RegionExprNode() {
    // children are contained in RegionExpr::all_region_expr_nodes,
    // and are deleted when ReleaseEvent is destroyed
  }

  bool has_binary_op() const {
    return
        op == RegionExprOperator::UNION ||
        op == RegionExprOperator::INTERSECT ||
        op == RegionExprOperator::DIFFERENCE;
  }

  RegionExprOperator op;

  region_id_t region_id;
  geometry_object_id_t geometry_object_id;

  RegionExprNode* left;
  RegionExprNode* right;

  void dump(const World* world = nullptr) const; // does not print any newlines
  std::string to_string(const World* world = nullptr, const bool for_datamodel = false) const;
};


// constructor and container for all region expr nodes
class RegionExpr {
public:
  RegionExpr() :
    root(nullptr) {
  }
  RegionExpr(const RegionExpr& other) {
    *this = other;
  }

  RegionExpr& operator=(const RegionExpr& other);

  ~RegionExpr();

  RegionExprNode* create_new_expr_node_leaf_surface_region(const region_id_t id);
  RegionExprNode* create_new_expr_node_leaf_geometry_object(const geometry_object_id_t id);

  RegionExprNode* create_new_region_expr_node_op(const RegionExprOperator op, RegionExprNode* left, RegionExprNode* right);

  // used when release_shape is ReleaseShape::Region
  RegionExprNode* root;

private:
  std::vector<RegionExprNode*> all_region_expr_nodes;
};

} /* namespace MCell */

#endif /* SRC4_REGION_EXPR_H_ */
