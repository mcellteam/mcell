/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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

std::shared_ptr<Observables> GenObservables::copy_observables() const {
  std::shared_ptr<Observables> res = std::make_shared<Observables>(DefaultCtorArgType());
  res->viz_outputs = viz_outputs;
  res->counts = counts;

  return res;
}

std::shared_ptr<Observables> GenObservables::deepcopy_observables(py::dict) const {
  std::shared_ptr<Observables> res = std::make_shared<Observables>(DefaultCtorArgType());
  for (const auto& item: viz_outputs) {
    res->viz_outputs.push_back((is_set(item)) ? item->deepcopy_viz_output() : nullptr);
  }
  for (const auto& item: counts) {
    res->counts.push_back((is_set(item)) ? item->deepcopy_count() : nullptr);
  }

  return res;
}

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

std::string GenObservables::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << "Observables" << ": " <<
      "\n" << ind + "  " << "viz_outputs=" << vec_ptr_to_str(viz_outputs, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "counts=" << vec_ptr_to_str(counts, all_details, ind + "  ");
  return ss.str();
}

py::class_<Observables> define_pybinding_Observables(py::module& m) {
  return py::class_<Observables, std::shared_ptr<Observables>>(m, "Observables", "Container used to hold observables-related model data. \nObservables are the measured values of the system. \nThis class also includes information on visualization of simulation.\n")
      .def(
          py::init<
            const std::vector<std::shared_ptr<VizOutput>>,
            const std::vector<std::shared_ptr<Count>>
          >(),
          py::arg("viz_outputs") = std::vector<std::shared_ptr<VizOutput>>(),
          py::arg("counts") = std::vector<std::shared_ptr<Count>>()
      )
      .def("__copy__", &Observables::copy_observables)
      .def("__deepcopy__", &Observables::deepcopy_observables, py::arg("memo"))
      .def("__str__", &Observables::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Observables::__eq__, py::arg("other"))
      .def("add_viz_output", &Observables::add_viz_output, py::arg("viz_output"), "Adds a reference to the viz_output object to the list of visualization output specifications.\n- viz_output\n")
      .def("add_count", &Observables::add_count, py::arg("count"), "Adds a reference to the count object to the list of count specifications.\n- count\n")
      .def("find_count", &Observables::find_count, py::arg("name"), "Finds a count object by its name, returns None if no such count is present.\n- name\n")
      .def("load_bngl_observables", &Observables::load_bngl_observables, py::arg("file_name"), py::arg("observables_path_or_file") = STR_UNSET, py::arg("parameter_overrides") = std::map<std::string, double>(), py::arg("observables_output_format") = CountOutputFormat::AUTOMATIC_FROM_EXTENSION, "Loads section observables from a BNGL file and creates Count objects according to it.\nAll elementary molecule types used in the seed species section must be defined in subsystem.\n\n- file_name: Path to the BNGL file.\n\n- observables_path_or_file: Directory prefix or file name where observable values will be stored.\nIf a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty \nstring/unset is used, each observable gets its own file and the output file format for created Count \nobjects is CountOutputFormat.DAT.\nWhen not set, this path is used: './react_data/seed_' + str(model.config.seed).zfill(5) + '/'.\nIf a file has a .gdat extension such as \n'./react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat', all observable are stored in this \nfile and the output file format for created Count objects is CountOutputFormat.GDAT.\nMust not be empty when observables_output_format is explicitly set to CountOutputFormat.GDAT.\n\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n- observables_output_format: Selection of output format. Default setting uses automatic detection\nbased on contents of the 'observables_path_or_file' attribute.\n             \n\n\n")
      .def("dump", &Observables::dump)
      .def_property("viz_outputs", &Observables::get_viz_outputs, &Observables::set_viz_outputs, py::return_value_policy::reference, "List of visualization outputs to be included in the model.\nThere is usually just one VizOutput object.   \n")
      .def_property("counts", &Observables::get_counts, &Observables::set_counts, py::return_value_policy::reference, "List of counts to be included in the model.\n")
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

