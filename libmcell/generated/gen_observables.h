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

#ifndef API_GEN_OBSERVABLES_H
#define API_GEN_OBSERVABLES_H

#include "api/common.h"

namespace MCell {
namespace API {

class Observables;
class Count;
class Subsystem;
class VizOutput;

class GenObservables {
public:
  virtual ~GenObservables() {}
  virtual bool __eq__(const Observables& other) const;
  virtual bool eq_nonarray_attributes(const Observables& other, const bool ignore_name = false) const;
  bool operator == (const Observables& other) const { return __eq__(other);}
  bool operator != (const Observables& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const ;

  // --- attributes ---
  std::vector<std::shared_ptr<VizOutput>> viz_outputs;
  virtual void set_viz_outputs(const std::vector<std::shared_ptr<VizOutput>> new_viz_outputs_) {
    viz_outputs = new_viz_outputs_;
  }
  virtual std::vector<std::shared_ptr<VizOutput>> get_viz_outputs() const {
    return viz_outputs;
  }

  std::vector<std::shared_ptr<Count>> counts;
  virtual void set_counts(const std::vector<std::shared_ptr<Count>> new_counts_) {
    counts = new_counts_;
  }
  virtual std::vector<std::shared_ptr<Count>> get_counts() const {
    return counts;
  }

  // --- methods ---
  virtual void add_viz_output(std::shared_ptr<VizOutput> viz_output) = 0;
  virtual void add_count(std::shared_ptr<Count> count) = 0;
  virtual std::shared_ptr<Count> find_count(const std::string& name) = 0;
  virtual void load_bngl_observables(const std::string& file_name, std::shared_ptr<Subsystem> subsystem, const std::string& output_files_prefix = "", const std::map<std::string, float_t>& parameter_overrides = std::map<std::string, float_t>()) = 0;
}; // GenObservables

class Observables;
py::class_<Observables> define_pybinding_Observables(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_OBSERVABLES_H
