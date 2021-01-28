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

#include "api/api_utils.h"

#include <time.h>
#include <string>

using namespace std;

namespace MCell {
namespace API {

Orientation convert_mcell_orientation(const orientation_t o) {
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


orientation_t convert_api_orientation(const Orientation o, const bool allow_default, const bool is_vol) {
  switch (o) {
    case Orientation::DEFAULT:
      if (!allow_default) {
        throw ValueError("Invalid Orientation value " + to_string((int)o) + " (DEFAULT).");
      }
      if (is_vol) {
        return ORIENTATION_NONE;
      }
      else {
        return ORIENTATION_UP;
      }
    case Orientation::DOWN:
      return ORIENTATION_DOWN;
    case Orientation::NONE:
      return ORIENTATION_NONE;
    case Orientation::UP:
      return ORIENTATION_UP;
    case Orientation::NOT_SET:
      throw ValueError("Invalid Orientation value " + to_string((int)o) + ".");
    case Orientation::ANY:
      return ORIENTATION_NONE;
    default:
      throw ValueError("Invalid Orientation value " + to_string((int)o) + ".");
  }
}


bool is_simple_species(const std::string& name) {
  // complex species always contain a parenthesis in their name
  return name.find('(') == std::string::npos;
}


// get current date/time, format is YYYY-MM-DD HH:mm:ss
const std::string get_current_date_time() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}

} // namespace API
} // namespace MCell

