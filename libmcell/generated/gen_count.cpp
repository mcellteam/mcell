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
#include "gen_count.h"
#include "api/count.h"
#include "api/count_term.h"

namespace MCell {
namespace API {

void GenCount::check_semantics() const {
}

void GenCount::set_initialized() {
  if (is_set(expression)) {
    expression->set_initialized();
  }
  initialized = true;
}

void GenCount::set_all_attributes_as_default_or_unset() {
  class_name = "Count";
  name = STR_UNSET;
  file_name = STR_UNSET;
  expression = nullptr;
  multiplier = 1;
  every_n_timesteps = 1;
  output_format = CountOutputFormat::AUTOMATIC_FROM_EXTENSION;
}

std::shared_ptr<Count> GenCount::copy_count() const {
  std::shared_ptr<Count> res = std::make_shared<Count>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->file_name = file_name;
  res->expression = expression;
  res->multiplier = multiplier;
  res->every_n_timesteps = every_n_timesteps;
  res->output_format = output_format;

  return res;
}

std::shared_ptr<Count> GenCount::deepcopy_count(py::dict) const {
  std::shared_ptr<Count> res = std::make_shared<Count>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->file_name = file_name;
  res->expression = is_set(expression) ? expression->deepcopy_count_term() : nullptr;
  res->multiplier = multiplier;
  res->every_n_timesteps = every_n_timesteps;
  res->output_format = output_format;

  return res;
}

bool GenCount::__eq__(const Count& other) const {
  return
    name == other.name &&
    file_name == other.file_name &&
    (
      (is_set(expression)) ?
        (is_set(other.expression) ?
          (expression->__eq__(*other.expression)) : 
          false
        ) :
        (is_set(other.expression) ?
          false :
          true
        )
     )  &&
    multiplier == other.multiplier &&
    every_n_timesteps == other.every_n_timesteps &&
    output_format == other.output_format;
}

bool GenCount::eq_nonarray_attributes(const Count& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    file_name == other.file_name &&
    (
      (is_set(expression)) ?
        (is_set(other.expression) ?
          (expression->__eq__(*other.expression)) : 
          false
        ) :
        (is_set(other.expression) ?
          false :
          true
        )
     )  &&
    multiplier == other.multiplier &&
    every_n_timesteps == other.every_n_timesteps &&
    output_format == other.output_format;
}

std::string GenCount::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "file_name=" << file_name << ", " <<
      "\n" << ind + "  " << "expression=" << "(" << ((expression != nullptr) ? expression->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "multiplier=" << multiplier << ", " <<
      "every_n_timesteps=" << every_n_timesteps << ", " <<
      "output_format=" << output_format;
  return ss.str();
}

py::class_<Count> define_pybinding_Count(py::module& m) {
  return py::class_<Count, std::shared_ptr<Count>>(m, "Count", "Represents a molecule or reaction count observable.\nWhat is counted is defined through a CounTerm tree and a reference to \nthe root of this tree is stored into attribute expression. \nThis tree usually contains just one node. \n  \n")
      .def(
          py::init<
            const std::string&,
            const std::string&,
            std::shared_ptr<CountTerm>,
            const double,
            const double,
            const CountOutputFormat
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("file_name") = STR_UNSET,
          py::arg("expression") = nullptr,
          py::arg("multiplier") = 1,
          py::arg("every_n_timesteps") = 1,
          py::arg("output_format") = CountOutputFormat::AUTOMATIC_FROM_EXTENSION
      )
      .def("check_semantics", &Count::check_semantics)
      .def("__copy__", &Count::copy_count)
      .def("__deepcopy__", &Count::deepcopy_count, py::arg("memo"))
      .def("__str__", &Count::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Count::__eq__, py::arg("other"))
      .def("get_current_value", &Count::get_current_value, "Returns the current value for this count. Cannot be used to count reactions.\nThe model must be initialized with this Count present as one of the observables.\n")
      .def("dump", &Count::dump)
      .def_property("name", &Count::get_name, &Count::set_name, "Name of a count may be specified when one needs to search for them later. \nWhen the count is created when a BNGL file is loaded, its name is set, for instance\nwhen the following BNGL code is loaded:\n\nbegin observables\n   Molecules Acount A\nend observables\n\nthe name is set to Acount.\n")
      .def_property("file_name", &Count::get_file_name, &Count::set_file_name, "File name where this observable values will be stored.\nFile extension or setting explicit output_format determines the output format.\nA) When not set, the value is set using seed during model initialization as follows: \nfile_name = './react_data/seed_' + str(model.config.seed).zfill(5) + '/' + name + '.dat'\nand the output format is set to CountOutputFormat.DAT in the constructor.\nB) When the file_name is set explicitly by the user and the extension is .dat such as here:\nfile_name = './react_data/seed_' + str(SEED).zfill(5) + '/' + name + '.dat'\nand the output format is set to CountOutputFormat.DAT in the constructor.\nFile names for individual Counts must be different.\nC) When the file_name is set explicitly by the user and the extension is .gdat such as here:\nfile_name = './react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat'\nand the output format is set to CountOutputFormat.GDAT in the constructor.\nThe file name is usually the same for all counts but one can \ncreate multiple gdat files with different observables.\nAll observables that are stored into a single .gdat file must have the same \nperiodicity specified by attribute every_n_timesteps.\nMust be set.\n")
      .def_property("expression", &Count::get_expression, &Count::set_expression, "The expression must be set to a root of an expression tree composed of CountTerms. \nIn the usual cases, there is just one CountTerm in this expression tree and its \nnode_type is ExprNodeType.LEAF.\nThe count expression tree defines CountTerm objects that are added or subtracted\nfrom each other.\n")
      .def_property("multiplier", &Count::get_multiplier, &Count::set_multiplier, "In some cases it might be useful to multiply the whole count by a constant to get \nfor instance concentration. The expression tree allows only addition and subtraction \nof count terms so such multiplication can be done through this attribute.\nIt can be also used to divide the resulting count by passing an inverse of the divisor (1/d).   \n")
      .def_property("every_n_timesteps", &Count::get_every_n_timesteps, &Count::set_every_n_timesteps, "Specifies periodicity of this count's output.\nValue is truncated (floored) to an integer.\nIf value is set to 0, this Count is used only on-demand through calls to its\nget_current_value method.  \n")
      .def_property("output_format", &Count::get_output_format, &Count::set_output_format, "Listed as the last attribute because the automatic default value\nis sufficient in most cases. \nSelection of output format. Default setting uses file extension  \nfrom attribute file_name. \nWhen set to CountOutputFormat.AUTOMATIC_FROM_EXTENSION, \nthis output_format is set automatically only in the Count's constructor. \n")
    ;
}

std::string GenCount::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("count") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("count")));
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
  ss << "m.Count(" << nl;
  if (name != STR_UNSET) {
    ss << ind << "name = " << "'" << name << "'" << "," << nl;
  }
  if (file_name != STR_UNSET) {
    ss << ind << "file_name = " << "'" << file_name << "'" << "," << nl;
  }
  if (is_set(expression)) {
    ss << ind << "expression = " << expression->export_to_python(out, ctx) << "," << nl;
  }
  if (multiplier != 1) {
    ss << ind << "multiplier = " << f_to_str(multiplier) << "," << nl;
  }
  if (every_n_timesteps != 1) {
    ss << ind << "every_n_timesteps = " << f_to_str(every_n_timesteps) << "," << nl;
  }
  if (output_format != CountOutputFormat::AUTOMATIC_FROM_EXTENSION) {
    ss << ind << "output_format = " << output_format << "," << nl;
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

} // namespace API
} // namespace MCell

