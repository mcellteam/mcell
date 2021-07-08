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

#ifndef LIBMCELL_API_BASE_EXPORT_CLASS_H_
#define LIBMCELL_API_BASE_EXPORT_CLASS_H_

#include "api/api_common.h"

namespace MCell {
namespace API {

class PythonExportContext;

// base class for all classes that provide export to python
class BaseExportClass {
public:
  virtual ~BaseExportClass() {
  }

  // may be overridden by the final class if needed
  virtual void set_all_custom_attributes_to_default() {
  }

  // used in generated export_to_python to optionally skip export of some objects
  virtual bool skip_python_export() const {
    return false;
  }

  // do not export vectors (used in Complex class)
  virtual bool skip_vectors_export() const {
    return false;
  }

  // export as a string without newlines (used in Complex class),
  // may not be applied to exported vectors
  virtual bool export_as_string_without_newlines() const {
    return false;
  }

  // do not check if already exported
  virtual bool export_even_if_already_exported() const {
    return false;
  }

  // either returns the whole definition as a string or prints definition to out and
  // returns name
  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx) {
    assert(false);
    return "Export to Python for a derived class is not implemented.";
  }
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_BASE_EXPORT_CLASS_H_ */
