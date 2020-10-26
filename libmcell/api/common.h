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

#ifndef LIBMCELL_API_COMMON_H
#define LIBMCELL_API_COMMON_H

#include <ostream>
#include <sstream>
#include <exception>
#include <string>
#include <memory>
#include <functional>
#include "pybind11/include/pybind11/pybind11.h" // make sure we won't include the system header
#include "pybind11/include/pybind11/functional.h"

namespace py = pybind11;
#include <vector>

#include "generated/gen_constants.h"
#include "generated/gen_names.h"

#include "defines.h"

namespace MCell {
namespace API {

// auxiliary method to simply convert to std::string for when concatenating string
static std::string S(const char* s) {
  return std::string(s);
}

const float_t FLT_UNSET = FLT_MAX;
const int INT_UNSET = INT32_MAX;
const long long LONG_UNSET = INT64_MAX;
const Vec3 VEC3_UNSET(POS_INVALID);
const void* const PTR_UNSET = nullptr;
const std::string STR_UNSET = "unset";

const std::string INTROSPECTED_OBJECT = "introspected object";

// specialization for each type
static inline bool is_set(const float_t a) {
  return a != FLT_UNSET;
}
static inline bool is_set(const int a) {
  return a != INT_UNSET; // Python does not use unsigned integers
}
static inline bool is_set(const long long a) {
  return a != LONG_UNSET; // Python does not use unsigned integers
}
static inline bool is_set(const Vec3& a) {
  return a != VEC3_UNSET;
}
static inline bool is_set(const void* a) {
  return a != PTR_UNSET;
}
static inline bool is_set(const std::string& a) {
  return a != STR_UNSET && a != "";
}
static inline bool is_set(const Orientation& a) {
  return a != Orientation::NOT_SET;
}
static inline bool is_set(const SurfacePropertyType& a) {
  return a != SurfacePropertyType::UNSET;
}
template<typename T>
static inline bool is_set(const std::shared_ptr<T>& a) {
  return a.use_count() != 0;
}
template<typename T>
static inline bool is_set(const std::vector<T>& a) {
  return !a.empty();
}


template<typename T1, typename T2>
int get_num_set(const T1& v1, const T2& v2) {
  int res = 0;
  if (is_set(v1)) {
    res++;
  }
  if (is_set(v2)) {
    res++;
  }
  return res;
}


template<typename T1, typename T2, typename T3>
uint get_num_set(const T1& v1, const T2& v2, const T3& v3) {
  int res = get_num_set(v1, v2);
  if (is_set(v3)) {
    res++;
  }
  return res;
}


template<typename T1, typename T2, typename T3, typename T4>
uint get_num_set(const T1& v1, const T2& v2, const T3& v3, const T4& v4) {
  int res = get_num_set(v1, v2, v3);
  if (is_set(v4)) {
    res++;
  }
  return res;
}


template<typename T1, typename T2, typename T3, typename T4, typename T5>
uint get_num_set(const T1& v1, const T2& v2, const T3& v3, const T4& v4, const T5& v5) {
  int res = get_num_set(v1, v2, v3, v4);
  if (is_set(v5)) {
    res++;
  }
  return res;
}


// NOTE: for some reason a definition in defines.h is not found, although it works fine for Vec3
static inline std::ostream & operator<<(std::ostream &out, const IVec3& a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}


template<typename T>
static inline std::string vec_ptr_to_str(const std::vector<T>& arr, const std::string ind="") {
  std::stringstream ss;
  ss << "[\n";
  for (size_t i = 0; i < arr.size(); i++) {
    ss << ind << " " << i << ":(" << arr[i]->to_str(ind) << ")";
    if (i + 1 != arr.size()) {
      ss << ", ";
    }
    ss << "\n";
  }
  ss << ind << "]";
  return ss.str();
}


template<typename T>
static inline std::string vec_nonptr_to_str(const std::vector<T>& arr, const std::string ind="") {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < arr.size(); i++) {
    ss << i << ":" << arr[i];
    if (i + 1 != arr.size()) {
      ss << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

template<typename T>
static inline std::string vec_nonptr_to_str(const std::vector<std::vector<T>>& arr, const std::string ind="") {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < arr.size(); i++) {
    ss << i << ":" << vec_nonptr_to_str(arr[i], ind + "  ");
    if (i + 1 != arr.size()) {
      ss << ", ";
    }
  }
  ss << "]";
  return ss.str();
}


template<typename T>
static inline bool vec_ptr_eq(const std::vector<std::shared_ptr<T>>& vec1, const std::vector<std::shared_ptr<T>>& vec2) {
  if (vec1.size() != vec2.size()) {
    return false;
  }
  for (size_t i = 0; i < vec1.size(); i++) {
    if (!vec1[i]->__eq__(*vec2[i])) {
      return false;
    }
  }
  return true;
}


template<typename T>
static inline void vec_set_initialized(const std::vector<std::shared_ptr<T>>& vec) {
  for (const std::shared_ptr<T>& item: vec) {
    item->set_initialized();
  }
}


// find item with name in a vector, returns shared null if not found
template<typename T>
static inline std::shared_ptr<T> vec_find_by_name(
    const std::vector<std::shared_ptr<T>>& vec, const std::string& name) {

  for (auto& item: vec) {
    if (item->name == name) {
      return item;
    }
  }
  return std::shared_ptr<T>(nullptr);
}


// Raised when an operation or function receives an argument that has the
// right type but an inappropriate value, and the situation is not described
// by a more precise exception such as IndexError.
typedef std::invalid_argument ValueError; // using naming from Python

// Raised when an error is detected that doesn’t fall in any of the other categories.
// The associated value is a string indicating what precisely went wrong.
typedef std::logic_error RuntimeError; // e.g. not defined?


} // namespace API
} // namespace MCell

#endif // LIBMCELL_API_COMMON_H
