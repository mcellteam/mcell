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

#include "api/python_export_utils.h"

namespace MCell {
namespace API {

bool PythonExportContext::already_exported(const BaseDataClass* obj) const {
  return exported_objects.count(obj) != 0;
}


std::string PythonExportContext::get_exported_name(const BaseDataClass* obj) const {
  assert(already_exported(obj));
  return exported_objects.find(obj)->second;
}


void PythonExportContext::add_exported(const BaseDataClass* obj, const std::string& name) {
  assert(!already_exported(obj));
  exported_objects[obj] = name;
}


uint PythonExportContext::postinc_counter(const std::string& underscored_class_name) {
  uint res;
  auto it = counters.find(underscored_class_name);
  if (it != counters.end()) {
    res = it->second;
    it->second++;
  }
  else {
    res = 0;
    counters[underscored_class_name] = 1;
  }

  return res;
}


std::string fix_id(const std::string& str) {
  std::string res;
  for (char c: str) {
    if (c == '+') {
      res += "_plus_";
    }
    else if (c == '-') {
      res += "_minus_";
    }
    else if (c == '?') {
      res += "_anybond_";
    }
    else if (c == '!') {
      res += "_bond_";
    }
    else if (c == '(') {
      res += "_ps_";
    }
    else if (c == ')') {
      res += "_pe_";
    }
    else if (
        c == ' ' || c == '.' || c == '_' ||
        c == ',' || c == '~') {
      res += "_";
    }
    else if (isalnum(c)) {
      res += c;
    }
    // ignoring the rest of the characters
  }
  return res;
}

} // namespace API
} // namespace MCell


