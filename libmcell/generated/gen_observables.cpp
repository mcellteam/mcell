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

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_observables.h"
#include "api/observables.h"
#include "api/count.h"
#include "api/viz_output.h"

namespace MCell {
namespace API {

bool GenObservables::__eq__(const Observables& other) const {
  return
    vec_ptr_eq(viz_outputs, other.viz_outputs) &&
    vec_ptr_eq(counts, other.counts);
}

bool GenObservables::eq_nonarray_attributes(const Observables& other, const bool ignore_name) const {
  return
    true /*viz_outputs*/ &&
    true /*counts*/;
}

std::string GenObservables::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Observables" << ": " <<
      "\n" << ind + "  " << "viz_outputs=" << vec_ptr_to_str(viz_outputs, ind + "  ") << ", " << "\n" << ind + "  " <<
      "counts=" << vec_ptr_to_str(counts, ind + "  ");
  return ss.str();
}

py::class_<Observables> define_pybinding_Observables(py::module& m) {
  return py::class_<Observables, std::shared_ptr<Observables>>(m, "Observables", "Neither VizOutput, nor Count have name, therefore there are no find_* methods.\n")
      .def(
          py::init<
            const std::vector<std::shared_ptr<VizOutput>>,
            const std::vector<std::shared_ptr<Count>>
          >(),
          py::arg("viz_outputs") = std::vector<std::shared_ptr<VizOutput>>(),
          py::arg("counts") = std::vector<std::shared_ptr<Count>>()
      )
      .def("__str__", &Observables::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Observables::__eq__, py::arg("other"))
      .def("add_viz_output", &Observables::add_viz_output, py::arg("viz_output"))
      .def("add_count", &Observables::add_count, py::arg("count"))
      .def("find_count", &Observables::find_count, py::arg("name"))
      .def("load_bngl_observables", &Observables::load_bngl_observables, py::arg("file_name"), py::arg("output_files_prefix") = "", py::arg("parameter_overrides") = std::map<std::string, float_t>(), "Loads section observables from a BNGL file and creates Count objects according to it.\nAll elementary molecule types used in the seed species section must be defined in subsystem.\n")
      .def("dump", &Observables::dump)
      .def_property("viz_outputs", &Observables::get_viz_outputs, &Observables::set_viz_outputs)
      .def_property("counts", &Observables::get_counts, &Observables::set_counts)
    ;
}

std::string GenObservables::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  std::string exported_name = "observables";

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Observables(" << nl;
  if (viz_outputs != std::vector<std::shared_ptr<VizOutput>>() && !skip_vectors_export()) {
    ss << ind << "viz_outputs = " << export_vec_viz_outputs(out, ctx, exported_name) << "," << nl;
  }
  if (counts != std::vector<std::shared_ptr<Count>>() && !skip_vectors_export()) {
    ss << ind << "counts = " << export_vec_counts(out, ctx, exported_name) << "," << nl;
  }
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

std::string GenObservables::export_vec_viz_outputs(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_viz_outputs";
  }
  else {
    exported_name = "viz_outputs";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < viz_outputs.size(); i++) {
    const auto& item = viz_outputs[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenObservables::export_vec_counts(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_counts";
  }
  else {
    exported_name = "counts";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < counts.size(); i++) {
    const auto& item = counts[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

} // namespace API
} // namespace MCell

