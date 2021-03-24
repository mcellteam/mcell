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

#ifndef API_WARNINGS_H
#define API_WARNINGS_H

#include "generated/gen_warnings.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class Warnings: public GenWarnings {
public:
  WARNINGS_CTOR()

  // must be called also manually during model initialization
  void check_semantics() const override {
    if (high_reaction_probability == WarningLevel::ERROR) {
      throw ValueError(S(NAME_CLASS_WARNINGS) + "." + NAME_HIGH_REACTION_PROBABILITY + " must not be set to " +
          NAME_ENUM_WARNING_LEVEL + "." + NAME_EV_ERROR + ".");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_WARNINGS_H
