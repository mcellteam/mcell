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

#ifndef API_RELEASE_PATTERN_H
#define API_RELEASE_PATTERN_H

#include "generated/gen_release_pattern.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class ReleasePattern: public GenReleasePattern {
public:
  RELEASE_PATTERN_CTOR()

  void check_semantics() const override {
    GenReleasePattern::check_semantics();

    if (train_interval < train_duration) {
      throw ValueError(S(NAME_RELEASE_PATTERN) + name + ": " +
          NAME_TRAIN_INTERVAL + " is shorter than " + NAME_TRAIN_DURATION + ".");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_RELEASE_PATTERN_H
