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

#ifndef API_SPECIES_H
#define API_SPECIES_H

#include <string>

#include "../generated/gen_species.h"
#include "common.h"

namespace MCell {
namespace API {

class Species: public GenSpecies {
public:
  SPECIES_CTOR()

  // actual manual implementation of a semantic check
  SemRes check_semantics(std::ostream& out) const override {
    SemRes base_res = GenSpecies::check_semantics(out);
    if (base_res != SemRes::OK) {
      return base_res;
    }

    if (is_set(diffusion_constant_2d) && is_set(diffusion_constant_3d)) {
      out << get_object_name() << "Only either 'diffusion_constant_2d' or 'diffusion_constant_3d' can be set.";
      return SemRes::ERROR;
    }

    return SemRes::OK;
  }
};

} // namespace API
} // namespace MCell

#endif // API_SPECIES_H
