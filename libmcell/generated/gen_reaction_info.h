/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_GEN_REACTION_INFO_H
#define API_GEN_REACTION_INFO_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class ReactionInfo;
class GeometryObject;
class ReactionRule;
class PythonExportContext;

class GenReactionInfo {
public:
  GenReactionInfo() {
  }
  GenReactionInfo(DefaultCtorArgType) {
  }
  virtual ~GenReactionInfo() {}
  std::shared_ptr<ReactionInfo> copy_reaction_info() const;
  std::shared_ptr<ReactionInfo> deepcopy_reaction_info(py::dict = py::dict()) const;
  virtual bool __eq__(const ReactionInfo& other) const;
  virtual bool eq_nonarray_attributes(const ReactionInfo& other, const bool ignore_name = false) const;
  bool operator == (const ReactionInfo& other) const { return __eq__(other);}
  bool operator != (const ReactionInfo& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const ;

  // --- attributes ---
  ReactionType type;
  virtual void set_type(const ReactionType new_type_) {
    type = new_type_;
  }
  virtual ReactionType get_type() const {
    return type;
  }

  std::vector<int> reactant_ids;
  virtual void set_reactant_ids(const std::vector<int> new_reactant_ids_) {
    reactant_ids = new_reactant_ids_;
  }
  virtual std::vector<int>& get_reactant_ids() {
    return reactant_ids;
  }

  std::vector<int> product_ids;
  virtual void set_product_ids(const std::vector<int> new_product_ids_) {
    product_ids = new_product_ids_;
  }
  virtual std::vector<int>& get_product_ids() {
    return product_ids;
  }

  std::shared_ptr<ReactionRule> reaction_rule;
  virtual void set_reaction_rule(std::shared_ptr<ReactionRule> new_reaction_rule_) {
    reaction_rule = new_reaction_rule_;
  }
  virtual std::shared_ptr<ReactionRule> get_reaction_rule() const {
    return reaction_rule;
  }

  double time;
  virtual void set_time(const double new_time_) {
    time = new_time_;
  }
  virtual double get_time() const {
    return time;
  }

  std::vector<double> pos3d;
  virtual void set_pos3d(const std::vector<double> new_pos3d_) {
    pos3d = new_pos3d_;
  }
  virtual std::vector<double>& get_pos3d() {
    return pos3d;
  }

  std::shared_ptr<GeometryObject> geometry_object;
  virtual void set_geometry_object(std::shared_ptr<GeometryObject> new_geometry_object_) {
    geometry_object = new_geometry_object_;
  }
  virtual std::shared_ptr<GeometryObject> get_geometry_object() const {
    return geometry_object;
  }

  int wall_index;
  virtual void set_wall_index(const int new_wall_index_) {
    wall_index = new_wall_index_;
  }
  virtual int get_wall_index() const {
    return wall_index;
  }

  std::vector<double> pos2d;
  virtual void set_pos2d(const std::vector<double> new_pos2d_) {
    pos2d = new_pos2d_;
  }
  virtual std::vector<double>& get_pos2d() {
    return pos2d;
  }

  // --- methods ---
}; // GenReactionInfo

class ReactionInfo;
py::class_<ReactionInfo> define_pybinding_ReactionInfo(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_REACTION_INFO_H
