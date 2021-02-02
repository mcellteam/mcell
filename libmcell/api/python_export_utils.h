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

#include <exception>
#include <string>
#include <vector>

#include "defines.h"
#include "api/python_export_constants.h"
#include "json/json.h"

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

typedef std::invalid_argument ConversionError;

// replace all characters that cannot be present in an identifier
std::string fix_id(const std::string& str);

std::string get_filename(const std::string& output_files_prefix, const std::string file_suffix, const char* ext);

void open_and_check_file_w_prefix(
    const std::string& output_files_prefix, const std::string file_suffix, std::ofstream& out,
    const bool for_append = false, const bool bngl = false);

std::string make_id(const std::string& s);

std::string fix_param_id(const std::string& str);

static std::string make_section_comment(const std::string text) {
  return BLOCK_BEGIN1 + text + BLOCK_BEGIN2 + "\n";
}

static std::string make_start_block_comment(const std::string text) {
  return BLOCK_BEGIN1 + text + BLOCK_BEGIN2;
}

static std::string make_end_block_comment(const std::string text) {
  return BLOCK_END1 + text + BLOCK_END2 + "\n";
}

static std::string make_enum_value(const std::string enum_name, const std::string value) {
  return MDOT + enum_name + "." + value;
}

static void gen_comma(std::ostream& out, Json::Value::ArrayIndex index, Json::Value& array) {
  if (index + 1 != array.size()) {
    out << ", ";
  }
}


template<typename T>
static void gen_comma(std::ostream& out, size_t index, const std::vector<T>& array) {
  if (index + 1 != array.size()) {
    out << ", ";
  }
}


template<typename T>
void print_comma(std::ostream& out, size_t index, const std::vector<T>& array) {
  if (index + 1 != array.size()) {
    out << ", ";
  }
}



// name might be empty
void gen_ctor_call(std::ostream& out, std::string name, std::string class_name, bool has_params = true);


static void gen_method_call(std::ostream& out, std::string obj, std::string method, std::string param = "") {
  out << obj << "." << method << "(" << param << ")\n";
}


template<typename T>
static void gen_param(std::ostream& out, std::string name, T value, bool comma) {
  out << IND << name << " = " << value << (comma?",":"") << "\n";
}


template<>
void gen_param(std::ostream& out, std::string name, Json::Value& value, bool comma) {
  out << IND << name << " = '" << value.asString() << "'" << (comma?",":"") << "\n";
}


template<>
void gen_param(std::ostream& out, std::string name, std::string value, bool comma) {
  out << IND << name << " = '" << value << "'" << (comma?",":"") << "\n";
}


template<>
void gen_param(std::ostream& out, std::string name, const char* const value, bool comma) {
  out << IND << name << " = '" << value << "'" << (comma?",":"") << "\n";
}


template<>
void gen_param(std::ostream& out, std::string name, int value, bool comma) {
  out << IND << name << " = " << value << (comma?",":"") << "\n";
}


template<>
void gen_param(std::ostream& out, std::string name, double value, bool comma) {
  out << IND << name << " = " << value << (comma?",":"") << "\n";
}

template<>
void gen_param(std::ostream& out, std::string name, bool value, bool comma) {
  out << IND << name << " = " << (value ? "True" : "False") << (comma?",":"") << "\n";
}

static void gen_param_id(std::ostream& out, std::string name, std::string id, bool comma) {
  out << IND << name << " = " << make_id(id) << (comma?",":"") << "\n";
}


static void gen_param_id(std::ostream& out, std::string name, Json::Value& id, bool comma) {
  out << IND << name << " = " << make_id(id.asString()) << (comma?",":"") << "\n";
}


void gen_param_expr(std::ostream& out, std::string name, const std::string& value, bool comma);

// this should be used when printing out floating point values (doubles)
static void gen_param_expr(std::ostream& out, std::string name, Json::Value& value, bool comma) {
  gen_param_expr(out, name, value.asString(), comma);
}


static void gen_param_enum(std::ostream& out, std::string name, std::string enum_name, std::string enum_value, bool comma) {
  out << IND << name << " = " << make_enum_value(enum_name, enum_value) << (comma?",":"") << "\n";
}


static void gen_param_list(std::ostream& out, std::string name, const std::vector<std::string>& values, bool comma, bool as_strings = false) {
  std::string q = (as_strings) ? "'" : "";

  out << IND << name << " = [";
  for (size_t i = 0; i < values.size(); i++) {
    out << q << values[i] << q;
    gen_comma(out, i , values);
  }

  out << "]" << (comma?",":"") << "\n";
}

static void gen_param_list_3_floats(std::ostream& out, std::string name, Json::Value& x, Json::Value& y, Json::Value& z, bool comma) {
  out << IND <<
      name << " = " <<
      "(" << x.asString() << ", " << y.asString() << ", " << z.asString() << ")" << (comma?",":"") << "\n";
}

template<typename T>
static void gen_assign(std::ostream& out, std::string var_name, T value) {
  out << var_name << " = " << value << "\n";
}


template<typename T>
static void gen_assign(std::ostream& out, std::string obj_name, std::string field_name1, T value) {
  out << obj_name << "." << field_name1 << " = " << value << "\n";
}


template<>
void gen_assign(std::ostream& out, std::string obj_name, std::string field_name1, bool value) {
  out << obj_name << "." << field_name1 << " = " << (value ? "True" : "False") << "\n";
}


static void gen_assign_str(std::ostream& out, std::string obj_name, std::string field_name1, std::string value) {
  out << obj_name << "." << field_name1 << " = '" << value << "'\n";
}


template<typename T>
static void gen_assign(std::ostream& out, std::string obj_name, std::string field_name1, std::string field_name2, T value) {
  out << obj_name << "." << field_name1 << "." << field_name2 << " = " << value << "\n";
}


// for some reason the template above casts double to int..
template<>
void gen_assign(std::ostream& out, std::string obj_name, std::string field_name1, std::string field_name2, double value) {
  out << obj_name << "." << field_name1 << "." << field_name2 << " = " << value << "\n";
}


template<>
void gen_assign(std::ostream& out, std::string obj_name, std::string field_name1, std::string field_name2, float value) {
  out << obj_name << "." << field_name1 << "." << field_name2 << " = " << value << "\n";
}


template<>
void gen_assign(std::ostream& out, std::string obj_name, std::string field_name1, std::string field_name2, bool value) {
  out << obj_name << "." << field_name1 << "." << field_name2 << " = " << (value ? "True" : "False") << "\n";
}

static void gen_assign_vec3(
    std::ostream& out, std::string obj_name, std::string field_name1, std::string field_name2,
    double x, double y, double z) {
  out << obj_name << "." << field_name1 << "." << field_name2 << " = [" << x << ", " << y << ", " << z  << "]\n";
}


static void gen_assign_vec3(
    std::ostream& out, std::string obj_name, std::string field_name1, std::string field_name2,
    float x, float y, float z) {
  out << obj_name << "." << field_name1 << "." << field_name2 << " = [" << x << ", " << y << ", " << z  << "]\n";
}

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_PYTHON_EXPORT_UTILS_H_ */
