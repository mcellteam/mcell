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

#ifndef LIBMCELL_API_BASE_INTROSPECTION_CLASS_H_
#define LIBMCELL_API_BASE_INTROSPECTION_CLASS_H_

#include "common.h"
#include "base_data_class.h"

namespace MCell {

class World;

namespace API {


// base class for all classes that hold the model input data
class BaseIntrospectionClass: public BaseDataClass {
public:
  BaseIntrospectionClass()
    : world(nullptr) {
    name = INTROSPECTED_OBJECT;
  }
  virtual ~BaseIntrospectionClass() {
  }

  void check_initialization() const {
    if (!initialized || world == nullptr) {
      throw RuntimeError(
          "Object of class " + class_name + " was not correctly initialized. "
          "Introspection objects canno be created independently. they must always be retrieved through "
          "methods of the " + NAME_CLASS_MODEL + " class."
      );
    }
  }

  void set_all_attributes_as_default_or_unset() {
    BaseDataClass::set_all_attributes_as_default_or_unset();
    world = nullptr;
  }

  // internal World pointer
  World* world;
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_BASE_INTROSPECTION_CLASS_H_ */
