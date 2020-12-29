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
#include "gen_chkpt_vol_mol.h"
#include "api/chkpt_vol_mol.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenChkptVolMol::check_semantics() const {
  if (!is_set(pos)) {
    throw ValueError("Parameter 'pos' must be set.");
  }
  if (!is_set(id)) {
    throw ValueError("Parameter 'id' must be set.");
  }
  if (!is_set(species)) {
    throw ValueError("Parameter 'species' must be set.");
  }
  if (!is_set(diffusion_time)) {
    throw ValueError("Parameter 'diffusion_time' must be set.");
  }
  if (!is_set(unimol_rx_time)) {
    throw ValueError("Parameter 'unimol_rx_time' must be set.");
  }
  if (!is_set(birthday)) {
    throw ValueError("Parameter 'birthday' must be set.");
  }
}

void GenChkptVolMol::set_initialized() {
  if (is_set(species)) {
    species->set_initialized();
  }
  initialized = true;
}

void GenChkptVolMol::set_all_attributes_as_default_or_unset() {
  class_name = "ChkptVolMol";
  pos = VEC3_UNSET;
  id = INT_UNSET;
  species = nullptr;
  diffusion_time = FLT_UNSET;
  unimol_rx_time = FLT_UNSET;
  birthday = FLT_UNSET;
}

bool GenChkptVolMol::__eq__(const ChkptVolMol& other) const {
  return
    pos == other.pos &&
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
    unimol_rx_time == other.unimol_rx_time &&
    birthday == other.birthday;
}

bool GenChkptVolMol::eq_nonarray_attributes(const ChkptVolMol& other, const bool ignore_name) const {
  return
    pos == other.pos &&
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
    unimol_rx_time == other.unimol_rx_time &&
    birthday == other.birthday;
}

std::string GenChkptVolMol::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "pos=" << pos << ", " <<
      "id=" << id << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "diffusion_time=" << diffusion_time << ", " <<
      "unimol_rx_time=" << unimol_rx_time << ", " <<
      "birthday=" << birthday;
  return ss.str();
}

py::class_<ChkptVolMol> define_pybinding_ChkptVolMol(py::module& m) {
  return py::class_<ChkptVolMol, BaseChkptMol, std::shared_ptr<ChkptVolMol>>(m, "ChkptVolMol")
      .def(
          py::init<
            const Vec3&,
            const int,
            std::shared_ptr<Species>,
            const float_t,
            const float_t,
            const float_t
          >(),
          py::arg("pos"),
          py::arg("id"),
          py::arg("species"),
          py::arg("diffusion_time"),
          py::arg("unimol_rx_time"),
          py::arg("birthday")
      )
      .def("check_semantics", &ChkptVolMol::check_semantics)
      .def("__str__", &ChkptVolMol::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ChkptVolMol::__eq__, py::arg("other"))
      .def("dump", &ChkptVolMol::dump)
      .def_property("pos", &ChkptVolMol::get_pos, &ChkptVolMol::set_pos)
    ;
}

std::string GenChkptVolMol::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "chkpt_vol_mol_" + std::to_string(ctx.postinc_counter("chkpt_vol_mol"));
  ctx.add_exported(this, exported_name);

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.ChkptVolMol(" << nl;
  ss << ind << "pos = " << "m.Vec3" << pos << "," << nl;
  ss << ind << "id = " << id << "," << nl;
  ss << ind << "species = " << species->export_to_python(out, ctx) << "," << nl;
  ss << ind << "diffusion_time = " << diffusion_time << "," << nl;
  ss << ind << "unimol_rx_time = " << unimol_rx_time << "," << nl;
  ss << ind << "birthday = " << birthday << "," << nl;
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

