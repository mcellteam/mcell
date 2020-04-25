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

#ifndef API_GEN_PARTITION_H
#define API_GEN_PARTITION_H

#include "../api/common.h"

namespace MCell {
namespace API {

#define PARTITION_CTOR() \
    Partition( \
        const float_t parition_dimension_ = 10, \
        const float_t subparition_dimension_ = 0.5 \
    ) { \
      class_name = "Partition"; \
      parition_dimension = parition_dimension_; \
      subparition_dimension = subparition_dimension_; \
      ctor_postprocess();\
    }

class GenPartition: public BaseDataClass {
public:
  void ctor_postprocess() override {}
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  float_t parition_dimension;
  virtual void set_parition_dimension(const float_t new_parition_dimension_) {
    parition_dimension = new_parition_dimension_;
  }
  virtual float_t get_parition_dimension() const {
    return parition_dimension;
  }

  float_t subparition_dimension;
  virtual void set_subparition_dimension(const float_t new_subparition_dimension_) {
    subparition_dimension = new_subparition_dimension_;
  }
  virtual float_t get_subparition_dimension() const {
    return subparition_dimension;
  }

  // --- methods ---
}; // GenPartition

class Partition;
py::class_<Partition> define_pybinding_Partition(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_PARTITION_H
