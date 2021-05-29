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
#include "gen_base_chkpt_mol.h"
#include "api/base_chkpt_mol.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenBaseChkptMol::check_semantics() const {
  if (!is_set(id)) {
    throw ValueError("Parameter 'id' must be set.");
  }
  if (!is_set(species)) {
    throw ValueError("Parameter 'species' must be set.");
  }
  if (!is_set(diffusion_time)) {
    throw ValueError("Parameter 'diffusion_time' must be set.");
  }
  if (!is_set(birthday)) {
    throw ValueError("Parameter 'birthday' must be set.");
  }
  if (!is_set(flags)) {
    throw ValueError("Parameter 'flags' must be set.");
  }
}

void GenBaseChkptMol::set_initialized() {
  if (is_set(species)) {
    species->set_initialized();
  }
  initialized = true;
}

void GenBaseChkptMol::set_all_attributes_as_default_or_unset() {
  class_name = "BaseChkptMol";
  id = INT_UNSET;
  species = nullptr;
  diffusion_time = FLT_UNSET;
  birthday = FLT_UNSET;
  flags = INT_UNSET;
  unimol_rx_time = FLT_UNSET;
}

std::shared_ptr<BaseChkptMol> GenBaseChkptMol::copy_base_chkpt_mol() const {
  if (initialized) {
    throw RuntimeError("Object of class BaseChkptMol cannot be cloned with 'copy' after this object was used in model initialization.");
  }

  std::shared_ptr<BaseChkptMol> res = std::make_shared<BaseChkptMol>(DefaultCtorArgType());
  res->class_name = class_name;
  res->id = id;
  res->species = species;
  res->diffusion_time = diffusion_time;
  res->birthday = birthday;
  res->flags = flags;
  res->unimol_rx_time = unimol_rx_time;

  return res;
}

std::shared_ptr<BaseChkptMol> GenBaseChkptMol::deepcopy_base_chkpt_mol(py::dict) const {
  if (initialized) {
    throw RuntimeError("Object of class BaseChkptMol cannot be cloned with 'deepcopy' after this object was used in model initialization.");
  }

  std::shared_ptr<BaseChkptMol> res = std::make_shared<BaseChkptMol>(DefaultCtorArgType());
  res->class_name = class_name;
  res->id = id;
  res->species = is_set(species) ? species->deepcopy_species() : nullptr;
  res->diffusion_time = diffusion_time;
  res->birthday = birthday;
  res->flags = flags;
  res->unimol_rx_time = unimol_rx_time;

  return res;
}

bool GenBaseChkptMol::__eq__(const BaseChkptMol& other) const {
  return
    id == other.id &&
    (
      (is_set(species)) ?
        (is_set(other.species) ?
          (species->__eq__(*other.species)) : 
          false
        ) :
        (is_set(other.species) ?
          false :
          true
        )
     )  &&
    diffusion_time == other.diffusion_time &&
    birthday == other.birthday &&
    flags == other.flags &&
    unimol_rx_time == other.unimol_rx_time;
}

bool GenBaseChkptMol::eq_nonarray_attributes(const BaseChkptMol& other, const bool ignore_name) const {
  return
    id == other.id &&
    (
      (is_set(species)) ?
        (is_set(other.species) ?
          (species->__eq__(*other.species)) : 
          false
        ) :
        (is_set(other.species) ?
          false :
          true
        )
     )  &&
    diffusion_time == other.diffusion_time &&
    birthday == other.birthday &&
    flags == other.flags &&
    unimol_rx_time == other.unimol_rx_time;
}

std::string GenBaseChkptMol::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "id=" << id << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "diffusion_time=" << diffusion_time << ", " <<
      "birthday=" << birthday << ", " <<
      "flags=" << flags << ", " <<
      "unimol_rx_time=" << unimol_rx_time;
  return ss.str();
}

py::class_<BaseChkptMol> define_pybinding_BaseChkptMol(py::module& m) {
  return py::class_<BaseChkptMol, std::shared_ptr<BaseChkptMol>>(m, "BaseChkptMol", "Base class for checkpointed molecules.\nNot to be used directly. All times are in seconds.\n")
      .def(
          py::init<
            const int,
            std::shared_ptr<Species>,
            const double,
            const double,
            const int,
            const double
          >(),
          py::arg("id"),
          py::arg("species"),
          py::arg("diffusion_time"),
          py::arg("birthday"),
          py::arg("flags"),
          py::arg("unimol_rx_time") = FLT_UNSET
      )
      .def("check_semantics", &BaseChkptMol::check_semantics)
      .def("__copy__", &BaseChkptMol::copy_base_chkpt_mol)
      .def("__deepcopy__", &BaseChkptMol::deepcopy_base_chkpt_mol, py::arg("memo"))
      .def("__str__", &BaseChkptMol::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &BaseChkptMol::__eq__, py::arg("other"))
      .def("dump", &BaseChkptMol::dump)
      .def_property("id", &BaseChkptMol::get_id, &BaseChkptMol::set_id)
      .def_property("species", &BaseChkptMol::get_species, &BaseChkptMol::set_species)
      .def_property("diffusion_time", &BaseChkptMol::get_diffusion_time, &BaseChkptMol::set_diffusion_time)
      .def_property("birthday", &BaseChkptMol::get_birthday, &BaseChkptMol::set_birthday)
      .def_property("flags", &BaseChkptMol::get_flags, &BaseChkptMol::set_flags)
      .def_property("unimol_rx_time", &BaseChkptMol::get_unimol_rx_time, &BaseChkptMol::set_unimol_rx_time)
    ;
}

std::string GenBaseChkptMol::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "base_chkpt_mol_" + std::to_string(ctx.postinc_counter("base_chkpt_mol"));
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
  ss << "m.BaseChkptMol(" << nl;
  ss << ind << "id = " << id << "," << nl;
  ss << ind << "species = " << species->export_to_python(out, ctx) << "," << nl;
  ss << ind << "diffusion_time = " << f_to_str(diffusion_time) << "," << nl;
  ss << ind << "birthday = " << f_to_str(birthday) << "," << nl;
  ss << ind << "flags = " << flags << "," << nl;
  if (unimol_rx_time != FLT_UNSET) {
    ss << ind << "unimol_rx_time = " << f_to_str(unimol_rx_time) << "," << nl;
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

