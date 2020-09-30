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

#ifndef UTILS_DATA_MODEL_TO_PYMCELL_BNG_GENERATOR_H_
#define UTILS_DATA_MODEL_TO_PYMCELL_BNG_GENERATOR_H_

#include <string>
#include <iostream>

#include "json/json.h"

namespace MCell {

class BNGGenerator {
public:
  BNGGenerator(std::ostream& bng_out_, Json::Value& mcell_)
    : bng_out(bng_out_), mcell(mcell_) {
  }

  void generate_parameters(std::ostream& python_out);

private:
  std::ostream& bng_out;
  Json::Value& mcell;
};

} /* namespace MCell */

#endif /* UTILS_DATA_MODEL_TO_PYMCELL_BNG_GENERATOR_H_ */
