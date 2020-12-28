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

#ifndef LIBMCELL_API_BASE_EXPORT_CLASS_H_
#define LIBMCELL_API_BASE_EXPORT_CLASS_H_

#include "common.h"

namespace MCell {
namespace API {

class PythonExportContext;

// base class for all classes that provide export to python
class BaseExportClass {
public:
  virtual ~BaseExportClass() {
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

  // either returns the whole definition as a string or prints definition to out and
  // returns name
  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx) const {
    assert(false);
    return "Export to Python for a derived class is not implemented.";
  }
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_BASE_EXPORT_CLASS_H_ */
