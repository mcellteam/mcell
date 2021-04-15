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
  return py::class_<Complex, std::shared_ptr<Complex>>(m, "Complex", "This class represents a complex molecule composed of molecule instances.\nIt is either defined using a BNGL string or using a list of elementary molecule instances.\nOn top of that, orientation may be defined.\nThis class is most often by calling its constructor as m.Complex(bngl_string) in cases where a \nfully qualified instance (such as for molecule releases) or a pattern (in observable counts) is needed.  \nComparison operator __eq__ first converts complexes to their canonical representation and \nthen does comparison so for instance m.Complex('A(b!1).B(a!1)') == m.Complex('B(a!2).A(b!2)').\n")
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
      .def("to_bngl_str", &Complex::to_bngl_str, "Creates a string that corresponds to its BNGL representation including compartments.")
      .def("as_species", &Complex::as_species, "Returns a Species object based on this Complex. All species-specific \nattributes are set to their default values and 'name' is set to value returned by \n'to_bngl_str()'.\n")
      .def("dump", &Complex::dump)
      .def_property("name", &Complex::get_name, &Complex::set_name, "When set, this complex instance is initialized from a BNGL string passed as this argument, \nthe string is parsed and elementary_molecules and compartment are initialized.\nOnly one of name or elementary_molecules can be set. \n")
      .def_property("elementary_molecules", &Complex::get_elementary_molecules, &Complex::set_elementary_molecules, py::return_value_policy::reference, "Individual molecule instances contained in the complex.\nOnly one of name or elementary_molecules can be set.\n")
      .def_property("orientation", &Complex::get_orientation, &Complex::set_orientation, "Specifies orientation of a molecule. \nWhen Orientation.DEFAULT if kept then during model initialization is\n'orientation' set to Orientation.NONE for volume complexes and to \nOrientation.UP for surface complexes.\nIgnored by derived class Species.\n")
      .def_property("compartment_name", &Complex::get_compartment_name, &Complex::set_compartment_name, "Specifies compartment name of this Complex.\nOnly one of 'orientation' and 'compartment_name' can be set. \nCorresponds to BNGL specification of a compartment for the whole complex '@COMP:'.\nIf a 2D/surface compartment is specified, the complex must be a surface complex and \norientation is set to Orientation.UP.\nIf a 3D/volume compartment is specified, the complex must be a volume complex and\norientation is set to Orientation.NONE.\nSets compartment to all elementary molecules whose compartment is unset. Does not override \nspecific compartments of elementary molecules that were already set.\nIf this is a volume complex (all elementary molecules have their diffusion_constant_3d set), \nall compartments of elementary molecules must be the same volume compartment.\nIf this is a surface complex (at least one elementary molecule has its their diffusion_constant_2d \nset), all compartments of surface elementary molecules must be the same, and\nall compartments of volume elementary molecules must be from the two neighboring \nvolume compartments.\n")
    ;
}

std::string GenComplex::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("complex") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("complex")));
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
  ss << "m.Complex(" << nl;
  if (name != STR_UNSET) {
    ss << ind << "name = " << "'" << name << "'" << "," << nl;
  }
  if (elementary_molecules != std::vector<std::shared_ptr<ElementaryMolecule>>() && !skip_vectors_export()) {
    ss << ind << "elementary_molecules = " << export_vec_elementary_molecules(out, ctx, exported_name) << "," << nl;
  }
  if (orientation != Orientation::DEFAULT) {
    ss << ind << "orientation = " << orientation << "," << nl;
  }
  if (compartment_name != STR_UNSET) {
    ss << ind << "compartment_name = " << "'" << compartment_name << "'" << "," << nl;
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

std::string GenComplex::export_vec_elementary_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
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

