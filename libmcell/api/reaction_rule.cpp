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

#include "api/complex.h"
#include "api/reaction_rule.h"
#include "bng/bngl_names.h"

using namespace std;

namespace MCell {
namespace API {

std::string ReactionRule::to_bngl_str() {
  string res;
  for (size_t i = 0; i < reactants.size(); i++) {
    res += reactants[i]->to_bngl_str();
    if (i + 1 != reactants.size()) {
      res += " + ";
    }
  }

  if (is_set(rev_rate)) {
    res += " <-> ";
  }
  else {
    res += " -> ";
  }

  if (!products.empty()) {
    for (size_t i = 0; i < products.size(); i++) {
      res += products[i]->to_bngl_str();
      if (i + 1 != products.size()) {
        res += " + ";
      }
    }
  }
  else {
    res += BNG::COMPLEX_Null;
  }

  return res;
}


} // namespace API
} // namespace MCell