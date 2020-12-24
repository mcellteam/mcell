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
#include "libs/pybind11/include/pybind11/stl.h"
#include "api/python_export_utils.h"
#include "gen_complex.h"
#include "api/complex.h"
#include "api/elementary_molecule.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenComplex::check_semantics() const {
}

void GenComplex::set_initialized() {
  vec_set_initialized(elementary_molecules);
  initialized = true;
}

void GenComplex::set_all_attributes_as_default_or_unset() {
  class_name = "Complex";
  name = STR_UNSET;
  elementary_molecules = std::vector<std::shared_ptr<ElementaryMolecule>>();
  orientation = Orientation::DEFAULT;
  compartment_name = STR_UNSET;
}

bool GenComplex::__eq__(const Complex& other) const {
  return
    name == other.name &&
    vec_ptr_eq(elementary_molecules, other.elementary_molecules) &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

bool GenComplex::eq_nonarray_attributes(const Complex& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*elementary_molecules*/ &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

std::string GenComplex::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "elementary_molecules=" << vec_ptr_to_str(elementary_molecules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation << ", " <<
      "compartment_name=" << compartment_name;
  return ss.str();
}

py::class_<Complex> define_pybinding_Complex(py::module& m) {
  return py::class_<Complex, std::shared_ptr<Complex>>(m, "Complex")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<ElementaryMolecule>>,
            const Orientation,
            const std::string&
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("elementary_molecules") = std::vector<std::shared_ptr<ElementaryMolecule>>(),
          py::arg("orientation") = Orientation::DEFAULT,
          py::arg("compartment_name") = STR_UNSET
      )
      .def("check_semantics", &Complex::check_semantics)
      .def("__str__", &Complex::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Complex::__eq__, py::arg("other"))
      .def("to_bngl_str", &Complex::to_bngl_str)
      .def("as_species", &Complex::as_species)
      .def("dump", &Complex::dump)
      .def_property("name", &Complex::get_name, &Complex::set_name)
      .def_property("elementary_molecules", &Complex::get_elementary_molecules, &Complex::set_elementary_molecules)
      .def_property("orientation", &Complex::get_orientation, &Complex::set_orientation)
      .def_property("compartment_name", &Complex::get_compartment_name, &Complex::set_compartment_name)
    ;
}

std::string GenComplex::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "complex_" + fix_id(name);
  ctx.add_exported(this, exported_name);

  std::stringstream ss;
  ss << exported_name << " = m.Complex(\n";
  if (name != STR_UNSET) {
    ss << "  name = " << "'" << name << "'" << ",\n";
  }
  if (elementary_molecules != std::vector<std::shared_ptr<ElementaryMolecule>>()) {
    ss << "  elementary_molecules = " << export_vec_elementary_molecules(out, ctx, exported_name) << ",\n";
  }
  if (orientation != Orientation::DEFAULT) {
    ss << "  orientation = " << orientation << ",\n";
  }
  if (compartment_name != STR_UNSET) {
    ss << "  compartment_name = " << "'" << name << "'" << ",\n";
  }
  ss << ")\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenComplex::export_vec_elementary_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < elementary_molecules.size(); i++) {
    const auto& item = elementary_molecules[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    std::string name = item->export_to_python(out, ctx);
    ss << name << ", ";
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

