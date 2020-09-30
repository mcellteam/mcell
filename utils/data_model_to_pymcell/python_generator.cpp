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

#include "generator_utils.h"
#include "python_generator.h"

using namespace std;
using namespace MCell::API;

namespace MCell {

using Json::Value;

void PythonGenerator::generate_single_parameter(ostream& out, Value& parameter) {
  out << "# " << parameter[KEY_PAR_DESCRIPTION].asString() << "\n";
  out << parameter[KEY_PAR_NAME].asString() << " = " << parameter[KEY_PAR_EXPRESSION].asString();
  string units = parameter[KEY_PAR_UNITS].asString();
  if (units != "") {
    out << " # units: " << units;
  }
  out << "\n\n";
}


void PythonGenerator::generate_parameters(std::ostream& out) {
  Value& parameter_system = get_node(mcell, KEY_PARAMETER_SYSTEM);
  if (parameter_system.isMember(KEY_MODEL_PARAMETERS)) {
    Value& parameter_list = get_node(parameter_system, KEY_MODEL_PARAMETERS);
    for (Value::ArrayIndex i = 0; i < parameter_list.size(); i++) {
      generate_single_parameter(out, parameter_list[i]);
    }
  }
  out << "\n";
}

} /* namespace MCell */
