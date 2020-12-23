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

#include <fstream>

#include "api/python_export_utils.h"
#include "api/python_export_constants.h"

namespace MCell {
namespace API {

using namespace std;

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


string get_filename(const string& output_files_prefix, const string file_suffix, const char* ext) {
  if (output_files_prefix == "" || output_files_prefix.back() == '/' || output_files_prefix.back() == '\\') {
    return output_files_prefix + file_suffix + ext;
  }
  else {
    return output_files_prefix + "_" + file_suffix + ext;
  }
}


void open_and_check_file_w_prefix(
    const string& output_files_prefix, const string file_suffix, ofstream& out,
    const bool for_append, const bool bngl) {

  string file_name = get_filename(output_files_prefix, file_suffix, (bngl) ? BNGL_EXT : PY_EXT);

  if (for_append) {
    cout << "Appending to " + file_name + ".\n";
    out.open(file_name, ios::app);
  }
  else {
    cout << "Generating file " + file_name + ".\n";
    out.open(file_name);
  }
  out.precision(FLOAT_OUT_PRECISION);
  if (!out.is_open()) {
    throw ConversionError("Could not open file '" + file_name + "' for writing.");
  }
}

} // namespace API
} // namespace MCell


