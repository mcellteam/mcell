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
#include "gen_viz_output.h"
#include "api/viz_output.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenVizOutput::check_semantics() const {
  if (!is_set(output_files_prefix)) {
    throw ValueError("Parameter 'output_files_prefix' must be set.");
  }
}

void GenVizOutput::set_initialized() {
  vec_set_initialized(species_list);
  initialized = true;
}

void GenVizOutput::set_all_attributes_as_default_or_unset() {
  class_name = "VizOutput";
  output_files_prefix = STR_UNSET;
  species_list = std::vector<std::shared_ptr<Species>>();
  mode = VizMode::ASCII;
  every_n_timesteps = 1;
}

std::shared_ptr<VizOutput> GenVizOutput::copy_viz_output() const {
  std::shared_ptr<VizOutput> res = std::make_shared<VizOutput>(DefaultCtorArgType());
  res->class_name = class_name;
  res->output_files_prefix = output_files_prefix;
  res->species_list = species_list;
  res->mode = mode;
  res->every_n_timesteps = every_n_timesteps;

  return res;
}

std::shared_ptr<VizOutput> GenVizOutput::deepcopy_viz_output(py::dict) const {
  std::shared_ptr<VizOutput> res = std::make_shared<VizOutput>(DefaultCtorArgType());
  res->class_name = class_name;
  res->output_files_prefix = output_files_prefix;
  for (const auto& item: species_list) {
    res->species_list.push_back((is_set(item)) ? item->deepcopy_species() : nullptr);
  }
  res->mode = mode;
  res->every_n_timesteps = every_n_timesteps;

  return res;
}

bool GenVizOutput::__eq__(const VizOutput& other) const {
  return
    output_files_prefix == other.output_files_prefix &&
    vec_ptr_eq(species_list, other.species_list) &&
    mode == other.mode &&
    every_n_timesteps == other.every_n_timesteps;
}

bool GenVizOutput::eq_nonarray_attributes(const VizOutput& other, const bool ignore_name) const {
  return
    output_files_prefix == other.output_files_prefix &&
    true /*species_list*/ &&
    mode == other.mode &&
    every_n_timesteps == other.every_n_timesteps;
}

std::string GenVizOutput::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "output_files_prefix=" << output_files_prefix << ", " <<
      "\n" << ind + "  " << "species_list=" << vec_ptr_to_str(species_list, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "mode=" << mode << ", " <<
      "every_n_timesteps=" << every_n_timesteps;
  return ss.str();
}

py::class_<VizOutput> define_pybinding_VizOutput(py::module& m) {
  return py::class_<VizOutput, std::shared_ptr<VizOutput>>(m, "VizOutput", "Defines a visualization output with locations of molecules \nthat can be then loaded by CellBlender.\n")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<Species>>,
            const VizMode,
            const double
          >(),
          py::arg("output_files_prefix"),
          py::arg("species_list") = std::vector<std::shared_ptr<Species>>(),
          py::arg("mode") = VizMode::ASCII,
          py::arg("every_n_timesteps") = 1
      )
      .def("check_semantics", &VizOutput::check_semantics)
      .def("__copy__", &VizOutput::copy_viz_output)
      .def("__deepcopy__", &VizOutput::deepcopy_viz_output, py::arg("memo"))
      .def("__str__", &VizOutput::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &VizOutput::__eq__, py::arg("other"))
      .def("dump", &VizOutput::dump)
      .def_property("output_files_prefix", &VizOutput::get_output_files_prefix, &VizOutput::set_output_files_prefix, "Prefix for the viz output files, the prefix value is computed from the simulation seed: \noutput_files_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene'.\n")
      .def_property("species_list", &VizOutput::get_species_list, &VizOutput::set_species_list, py::return_value_policy::reference, "Specifies a list of species to be visualized, when empty, all_species will be generated.")
      .def_property("mode", &VizOutput::get_mode, &VizOutput::set_mode, "Specified the output format of the visualization files. \nVizMode.ASCII is a readable representation, VizMode.CELLBLENDER is a binary representation \nthat cannot be read using a text editor but is faster to generate. \n")
      .def_property("every_n_timesteps", &VizOutput::get_every_n_timesteps, &VizOutput::set_every_n_timesteps, "Specifies periodicity of visualization output.\nValue is truncated (floored) to an integer.\nValue 0 means that the viz output is ran only once at iteration 0. \n")
    ;
}

std::string GenVizOutput::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "viz_output_" + std::to_string(ctx.postinc_counter("viz_output"));
  if (!export_even_if_already_exported()) {
    ctx.add_exported(this, exported_name);
  }

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.VizOutput(" << nl;
  ss << ind << "output_files_prefix = " << "'" << output_files_prefix << "'" << "," << nl;
  if (species_list != std::vector<std::shared_ptr<Species>>() && !skip_vectors_export()) {
    ss << ind << "species_list = " << export_vec_species_list(out, ctx, exported_name) << "," << nl;
  }
  if (mode != VizMode::ASCII) {
    ss << ind << "mode = " << mode << "," << nl;
  }
  if (every_n_timesteps != 1) {
    ss << ind << "every_n_timesteps = " << f_to_str(every_n_timesteps) << "," << nl;
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

std::string GenVizOutput::export_vec_species_list(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < species_list.size(); i++) {
    const auto& item = species_list[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

