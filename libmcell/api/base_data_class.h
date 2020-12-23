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

#ifndef LIBMCELL_API_BASE_DATA_CLASS_H_
#define LIBMCELL_API_BASE_DATA_CLASS_H_

#include "common.h"

namespace MCell {
namespace API {

class PythonExportContext;

// base class for all classes that hold the model input data
class BaseDataClass {
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
  virtual std::string to_str(const std::string ind="") const {
    assert(false);
    return "String dump for a derived class is not implemented.";
  }

  // initialization for custom constructors
  virtual void set_all_attributes_as_default_or_unset() {
    name = STR_UNSET;
    initialized = false;
    cached_data_are_uptodate = false;
  };

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx) const {
    assert(false);
    return "Export to Python for a derived class is not implemented.";
  }

  // calls virtual method, usually no need to override
  virtual void dump() const {
    std::cout << to_str() << "\n";
  }
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_BASE_DATA_CLASS_H_ */
