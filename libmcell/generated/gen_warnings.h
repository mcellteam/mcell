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

#ifndef API_GEN_WARNINGS_H
#define API_GEN_WARNINGS_H

#include "api/common.h"
#include "api/base_data_class.h"

namespace MCell {
namespace API {

class Warnings;

#define WARNINGS_CTOR() \
    Warnings( \
    ) { \
      class_name = "Warnings"; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenWarnings: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const Warnings& other) const;
  virtual bool eq_nonarray_attributes(const Warnings& other, const bool ignore_name = false) const;
  bool operator == (const Warnings& other) const { return __eq__(other);}
  bool operator != (const Warnings& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  std::string export_to_python(std::ostream& out) const override;


  // --- attributes ---
  // --- methods ---
}; // GenWarnings

class Warnings;
py::class_<Warnings> define_pybinding_Warnings(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_WARNINGS_H
