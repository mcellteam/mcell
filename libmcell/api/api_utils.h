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

#ifndef LIBMCELL_API_API_UTILS_H_
#define LIBMCELL_API_API_UTILS_H_

#include <memory>
#include <vector>

#include "api/api_common.h"

// functions defined here are not used in generated files

namespace MCell {
namespace API {

bool is_simple_species(const std::string& name);
Orientation convert_mcell_orientation(const orientation_t o);
orientation_t convert_api_orientation(const Orientation o, const bool allow_default = false, const bool is_vol = true);

static inline std::vector<double> mult_vec(const std::vector<double>& in, const double mult) {
  std::vector<double> res(in.size());
  for (size_t i = 0; i < in.size(); i++) {
    res[i] = in[i] * mult;
  }
  return res;
}

template<class T>
void append_to_vec(
    std::vector<std::shared_ptr<T>>& dst,
    const std::shared_ptr<T>& item,
    const bool allow_same_name_different_contents = false) {

  if (!allow_same_name_different_contents) {

    // check if item with this name already exists
    for (std::shared_ptr<T>& existing: dst) {
      if (item->name == existing->name) {
        // must be identical
        if (!item->__eq__(*existing)) {
          throw ValueError(
              "Adding object of " + item->class_name + " with name '" + item->name +
              "' caused an error, object with the same name is already present but it is different, " +
              "\nexisting:\n" + existing->to_str() + "\nvs. new\n:" + item->to_str() + ".\n"
              "-- end of error message reporting existence of two objects with identical names --"
          );
        }
        else {
          std::cerr << "Warning: adding of " + item->class_name + " with name '" + item->name +
              "' is ignored, identical object is already present.\n";
          return;
        }
      }
    }
  }

  dst.push_back(item);
}


template<class T>
void append_to_vec_canonical_name(
    std::vector<std::shared_ptr<T>>& dst,
    const std::shared_ptr<T>& item) {

  // check if item with this name already exists
  for (std::shared_ptr<T>& existing: dst) {
    bool are_equal = item->__eq__(*existing);

    std::string name_in_msg = "object of " + item->class_name + " with" +
        (is_set(item->name) ? " name '" + item->name + "' and" : S("")) +
        " canonical name '"  + item->get_canonical_name() + "'";

    if (item->get_canonical_name() == existing->get_canonical_name()) {
      // must be identical
      if (!are_equal) {
        throw ValueError(
            "Adding " + name_in_msg + " caused an error, object with the same canonical name is already present but it is different, "
            "\nexisting:\n" + existing->to_str() + "\nvs. new:\n " + item->to_str() + ".\n"
            "-- end of error message reporting existence of two objects with identical names --"
        );
      }
      else {
        if (item->warn_if_adding_identical_object()) {
          std::cerr << "Warning: adding of " + name_in_msg + " is ignored, object with the same canonical name is already present.\n";
        }
        return;
      }

    }

    // also check, if name is set, that I cannot have different objects with the same name
    if (is_set(item->name) && item->name == existing->name && !are_equal) {
      throw ValueError(
          "Adding object of " + name_in_msg + " caused an error, object with the same name is already present but it is different."
      );
    }
  }

  dst.push_back(item);
}


template<class T>
void append_vec_to_vec(
    std::vector<std::shared_ptr<T>>& dst,
    const std::vector<std::shared_ptr<T>>& src,
    const bool allow_same_name_different_contents = false
) {

  for (const std::shared_ptr<T>& item: src) {
    append_to_vec(
        dst, item, allow_same_name_different_contents);
  }
}


template<class T>
void append_vec_to_vec_canonical_name(
    std::vector<std::shared_ptr<T>>& dst,
    const std::vector<std::shared_ptr<T>>& src
) {

  for (const std::shared_ptr<T>& item: src) {
    append_to_vec_canonical_name(dst, item);
  }
}

} // namespace API
} // namespace MCell

#endif // LIBMCELL_API_API_UTILS_H_
