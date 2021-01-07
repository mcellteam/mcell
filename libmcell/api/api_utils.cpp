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

#include "api_utils.h"

#include <iomanip>

namespace MCell {
namespace API {

Orientation convert_orientation(const orientation_t o) {
  switch (o) {
    case ORIENTATION_DOWN:
      return Orientation::DOWN;
    case ORIENTATION_NONE:
      return Orientation::NONE;
    case ORIENTATION_UP:
      return Orientation::UP;
    case ORIENTATION_NOT_SET:
      return Orientation::NOT_SET;
    default:
      assert(false);
      return Orientation::NOT_SET;
  }
}


bool is_simple_species(const std::string& name) {
  // complex species always contain a parenthesis in their name
  return name.find('(') == std::string::npos;
}


std::string get_seed_dir_name(const int seed) {
  std::stringstream seed_num;
  seed_num << std::setfill('0') << std::setw(DEFAULT_SEED_DIR_DIGITS) << seed;
  return DEFAULT_SEED_DIR_PREFIX + seed_num.str();
}

} // namespace API
} // namespace MCell

