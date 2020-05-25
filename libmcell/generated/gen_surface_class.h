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

#ifndef API_GEN_SURFACE_CLASS_H
#define API_GEN_SURFACE_CLASS_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class Species;

#define SURFACE_CLASS_CTOR() \
    SurfaceClass( \
        const std::string& name_, \
        std::shared_ptr<Species> reflective_ = nullptr, \
        std::shared_ptr<Species> transparent_ = nullptr, \
        std::shared_ptr<Species> absorptive_ = nullptr \
    ) { \
      class_name = "SurfaceClass"; \
      name = name_; \
      reflective = reflective_; \
      transparent = transparent_; \
      absorptive = absorptive_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenSurfaceClass: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;

  bool __eq__(const GenSurfaceClass& other) const;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Species> reflective;
  virtual void set_reflective(std::shared_ptr<Species> new_reflective_) {
    if (initialized) {
      throw RuntimeError("Value 'reflective' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    reflective = new_reflective_;
  }
  virtual std::shared_ptr<Species> get_reflective() const {
    return reflective;
  }

  std::shared_ptr<Species> transparent;
  virtual void set_transparent(std::shared_ptr<Species> new_transparent_) {
    if (initialized) {
      throw RuntimeError("Value 'transparent' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    transparent = new_transparent_;
  }
  virtual std::shared_ptr<Species> get_transparent() const {
    return transparent;
  }

  std::shared_ptr<Species> absorptive;
  virtual void set_absorptive(std::shared_ptr<Species> new_absorptive_) {
    if (initialized) {
      throw RuntimeError("Value 'absorptive' of object with name " + name + " (class " + class_name + ")"
                         "cannot be set after model was initialized.");
    }
    absorptive = new_absorptive_;
  }
  virtual std::shared_ptr<Species> get_absorptive() const {
    return absorptive;
  }

  // --- methods ---
}; // GenSurfaceClass

class SurfaceClass;
py::class_<SurfaceClass> define_pybinding_SurfaceClass(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SURFACE_CLASS_H
