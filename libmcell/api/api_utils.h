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

#ifndef LIBMCELL_API_API_UTILS_H_
#define LIBMCELL_API_API_UTILS_H_

#include <memory>
#include <vector>

#include "common.h"

namespace MCell {
namespace API {

template<class T>
void append_to_vector(
    std::vector<std::shared_ptr<T>>& dst,
    const std::shared_ptr<T>& item) {

  // check if item with this name already exists
  for (std::shared_ptr<T>& existing: dst) {
    if (item->name == existing->name) {
      // must be identical
      if (item->__eq__(*existing)) {
        throw ValueError(
            "Adding object of " + item->class_name +
            " caused an error, object with the same name is already present but it is different."
        );
      }
    }
  }
}


template<class T>
void append_vector_to_vector(
    std::vector<std::shared_ptr<T>>& dst,
    const std::vector<std::shared_ptr<T>>& src) {

  for (const std::shared_ptr<T>& item: src) {
    append_to_vector(dst, item);
  }
}

} /* namespace API */
} /* namespace MCell */

#endif /* LIBMCELL_API_API_UTILS_H_ */
