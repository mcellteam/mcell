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

#include "api/api_common.h"
#include "api/base_export_class.h"

namespace MCell {
namespace API {

class Observables;
class Count;
class VizOutput;
class PythonExportContext;

#define OBSERVABLES_CTOR() \
    Observables( \
        const std::vector<std::shared_ptr<VizOutput>> viz_outputs_ = std::vector<std::shared_ptr<VizOutput>>(), \
        const std::vector<std::shared_ptr<Count>> counts_ = std::vector<std::shared_ptr<Count>>() \
    ) { \
      viz_outputs = viz_outputs_; \
      counts = counts_; \
    } \
    Observables(DefaultCtorArgType){ \
    }

class GenObservables: public BaseExportClass {
public:
  GenObservables() {
  }
  GenObservables(DefaultCtorArgType) {
  }
  virtual ~GenObservables() {}
  std::shared_ptr<Observables> copy_observables() const;
  std::shared_ptr<Observables> deepcopy_observables(py::dict = py::dict()) const;
  virtual bool __eq__(const Observables& other) const;
  virtual bool eq_nonarray_attributes(const Observables& other, const bool ignore_name = false) const;
  bool operator == (const Observables& other) const { return __eq__(other);}
  bool operator != (const Observables& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const ;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);
  virtual std::string export_vec_viz_outputs(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_counts(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::shared_ptr<VizOutput>> viz_outputs;
  virtual void set_viz_outputs(const std::vector<std::shared_ptr<VizOutput>> new_viz_outputs_) {
    viz_outputs = new_viz_outputs_;
  }
  virtual std::vector<std::shared_ptr<VizOutput>>& get_viz_outputs() {
    return viz_outputs;
  }

  std::vector<std::shared_ptr<Count>> counts;
  virtual void set_counts(const std::vector<std::shared_ptr<Count>> new_counts_) {
    counts = new_counts_;
  }
  virtual std::vector<std::shared_ptr<Count>>& get_counts() {
    return counts;
  }

  // --- methods ---
  virtual void add_viz_output(std::shared_ptr<VizOutput> viz_output) = 0;
  virtual void add_count(std::shared_ptr<Count> count) = 0;
  virtual std::shared_ptr<Count> find_count(const std::string& name) = 0;
  virtual void load_bngl_observables(const std::string& file_name, const std::string& output_files_prefix = "", const std::map<std::string, double>& parameter_overrides = std::map<std::string, double>()) = 0;
}; // GenObservables

class Observables;
py::class_<Observables> define_pybinding_Observables(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_OBSERVABLES_H
