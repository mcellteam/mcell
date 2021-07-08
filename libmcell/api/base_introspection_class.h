/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef LIBMCELL_API_BASE_INTROSPECTION_CLASS_H_
#define LIBMCELL_API_BASE_INTROSPECTION_CLASS_H_

#include "api/api_common.h"
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
    // - introspected objects are assumed to be initialized because they are returned
    //   by API methods and this does not depend on model initialization
    // - this flag is checked in set_* methods and must be true to avoid ignoring writes to attributes
    initialized = true;
  }
  virtual ~BaseIntrospectionClass() {
  }

  void check_initialization() const {
    if (world == nullptr) {
      throw RuntimeError(
          "Object of class " + class_name + " was not correctly initialized. "
          "Introspection objects cannot be created independently. they must always be retrieved through "
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
