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

#ifndef LIBMCELL_API_PYTHON_EXPORT_UTILS_H_
#define LIBMCELL_API_PYTHON_EXPORT_UTILS_H_

#include "defines.h"

namespace MCell {
namespace API {

class BaseDataClass;


// class used to hold data when exporting objects to Python, e.g. for checkpointing
class PythonExportContext {
public:
  bool already_exported(const BaseDataClass* obj) const;
  std::string get_exported_name(const BaseDataClass* obj) const;
  void add_exported(const BaseDataClass* obj, const std::string& name);

  // returns current counter value (starts at 0) and increments counter for the current
  // class
  uint postinc_counter(const std::string& underscored_class_name);
private:
  std::map<const BaseDataClass*, std::string> exported_objects;

  std::map<std::string, uint> counters;
};


// replace all characters that cannot be present in an identifier
std::string fix_id(const std::string& str);

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_PYTHON_EXPORT_UTILS_H_ */
