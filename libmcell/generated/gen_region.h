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

#ifndef API_GEN_REGION_H
#define API_GEN_REGION_H

#include "api/api_common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Region;
class PythonExportContext;

#define REGION_CTOR() \
    Region( \
        const RegionNodeType node_type_ = RegionNodeType::UNSET, \
        std::shared_ptr<Region> left_node_ = nullptr, \
        std::shared_ptr<Region> right_node_ = nullptr \
    ) { \
      class_name = "Region"; \
      node_type = node_type_; \
      left_node = left_node_; \
      right_node = right_node_; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Region(DefaultCtorArgType) : \
      GenRegion(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
      set_all_custom_attributes_to_default(); \
    }

class GenRegion: public BaseDataClass {
public:
  GenRegion() {
  }
  GenRegion(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  std::shared_ptr<Region> copy_region() const;
  std::shared_ptr<Region> deepcopy_region(py::dict = py::dict()) const;
  virtual bool __eq__(const Region& other) const;
  virtual bool eq_nonarray_attributes(const Region& other, const bool ignore_name = false) const;
  bool operator == (const Region& other) const { return __eq__(other);}
  bool operator != (const Region& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) override;


  // --- attributes ---
  RegionNodeType node_type;
  virtual void set_node_type(const RegionNodeType new_node_type_) {
    if (initialized) {
      throw RuntimeError("Value 'node_type' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    node_type = new_node_type_;
  }
  virtual RegionNodeType get_node_type() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return node_type;
  }

  std::shared_ptr<Region> left_node;
  virtual void set_left_node(std::shared_ptr<Region> new_left_node_) {
    if (initialized) {
      throw RuntimeError("Value 'left_node' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    left_node = new_left_node_;
  }
  virtual std::shared_ptr<Region> get_left_node() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return left_node;
  }

  std::shared_ptr<Region> right_node;
  virtual void set_right_node(std::shared_ptr<Region> new_right_node_) {
    if (initialized) {
      throw RuntimeError("Value 'right_node' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    right_node = new_right_node_;
  }
  virtual std::shared_ptr<Region> get_right_node() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return right_node;
  }

  // --- methods ---
  virtual std::shared_ptr<Region> __add__(std::shared_ptr<Region> other) = 0;
  virtual std::shared_ptr<Region> __sub__(std::shared_ptr<Region> other) = 0;
  virtual std::shared_ptr<Region> __mul__(std::shared_ptr<Region> other) = 0;
}; // GenRegion

class Region;
py::class_<Region> define_pybinding_Region(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_REGION_H
