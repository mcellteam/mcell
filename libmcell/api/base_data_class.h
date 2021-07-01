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

#ifndef LIBMCELL_API_BASE_DATA_CLASS_H_
#define LIBMCELL_API_BASE_DATA_CLASS_H_

#include "api/api_common.h"
#include "api/base_export_class.h"

namespace MCell {
namespace API {

class PythonExportContext;

// base class for all classes that hold the model input data
class BaseDataClass: public BaseExportClass {
public:
  BaseDataClass()
    : class_name(STR_UNSET), name(STR_UNSET), initialized(false), cached_data_are_uptodate(false) {
  }
  virtual ~BaseDataClass() {
  }

  // we are storing class name for reporting
  std::string class_name;
  virtual void set_class_name(const std::string& class_name_) {
    class_name = class_name_;
  }
  virtual const std::string& get_class_name() const {
    return class_name;
  }

  // every object defined by the MCell API might have its name
  std::string name;
  virtual void set_name(const std::string& name_) {
    cached_data_are_uptodate = false;
    name = name_;
  }
  virtual const std::string& get_name() const {
    cached_data_are_uptodate = false; // might be modified in theory
    return name;
  }

  bool initialized;
  virtual void set_initialized() = 0;

  // attribute used when some caching is employed
  mutable bool cached_data_are_uptodate;

  // this method is used to identify this particular object in error messages
  virtual std::string get_object_name() const {
    return get_class_name() + " '" + get_name() + "'";
  }

  // empty implementation, to be overridden in actual derived classes
  virtual void postprocess_in_ctor() { };

  // empty implementation, to be overridden in actual derived classes
  virtual void check_semantics() const { };

  // empty implementation, to be overridden in actual derived classes
  virtual std::string to_str(const bool all_details=false, const std::string ind="") const {
    assert(false);
    return "String dump for a derived class is not implemented.";
  }

  // initialization for custom constructors
  virtual void set_all_attributes_as_default_or_unset() {
    name = STR_UNSET;
    initialized = false;
    cached_data_are_uptodate = false;
  };

  // calls virtual method, usually no need to override
  virtual void dump() const {
    std::cout << to_str() << "\n";
  }
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_BASE_DATA_CLASS_H_ */
