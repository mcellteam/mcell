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

#ifndef SRC4_PYMCELLCONVERTER_H_
#define SRC4_PYMCELLCONVERTER_H_

#include <string>
#include <fstream>

// use different namespace?
namespace MCell {

const char* const GEOMETRY_SUFFIX = "geometry";
const char* const PY_EXT = ".py";

const char* const MCELL_IMPORT = "import mcell as m\n\n";

class PymcellConverter {
  void convert(World* world_);

private:
  void convert_geometry();

  std::string file_prefix;
};

} /* namespace MCell */

#endif /* SRC4_PYMCELLCONVERTER_H_ */
