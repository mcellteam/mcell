/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef API_RELEASE_SITE_H
#define API_RELEASE_SITE_H

#include <string>

#include "../generated/gen_release_site.h"
#include "common.h"

namespace MCell {
namespace API {

class ReleaseSite: public GenReleaseSite {
public:
  RELEASE_SITE_CTOR()

  // actual manual implementation of a semantic check
  SemRes check_semantics(std::ostream& out) const override {
    SemRes base_res = GenReleaseSite::check_semantics(out);
    if (base_res != SemRes::OK) {
      return base_res;
    }

    if (is_set(site_diameter) && is_set(site_radius)) {
      out << "Only either 'site_diameter' or 'site_radius' can be set.\n";
      return SemRes::ERROR;
    }

    return SemRes::OK;
  }
};



} // namespace API
} // namespace MCell

#endif // API_RELEASE_SITE_H
