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
  if (!is_set(birthday)) {
    throw ValueError("Parameter 'birthday' must be set.");
  }
  if (!is_set(flags)) {
    throw ValueError("Parameter 'flags' must be set.");
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
  birthday = FLT_UNSET;
  flags = INT_UNSET;
  unimol_rxn_time = FLT_UNSET;
}

std::shared_ptr<ChkptVolMol> GenChkptVolMol::copy_chkpt_vol_mol() const {
  std::shared_ptr<ChkptVolMol> res = std::make_shared<ChkptVolMol>(DefaultCtorArgType());
  res->class_name = class_name;
  res->pos = pos;
  res->id = id;
  res->species = species;
  res->diffusion_time = diffusion_time;
  res->birthday = birthday;
  res->flags = flags;
  res->unimol_rxn_time = unimol_rxn_time;

  return res;
}

std::shared_ptr<ChkptVolMol> GenChkptVolMol::deepcopy_chkpt_vol_mol(py::dict) const {
  std::shared_ptr<ChkptVolMol> res = std::make_shared<ChkptVolMol>(DefaultCtorArgType());
  res->class_name = class_name;
  res->pos = pos;
  res->id = id;
  res->species = is_set(species) ? species->deepcopy_species() : nullptr;
  res->diffusion_time = diffusion_time;
  res->birthday = birthday;
  res->flags = flags;
  res->unimol_rxn_time = unimol_rxn_time;

  return res;
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
    birthday == other.birthday &&
    flags == other.flags &&
    unimol_rxn_time == other.unimol_rxn_time;
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
    birthday == other.birthday &&
    flags == other.flags &&
    unimol_rxn_time == other.unimol_rxn_time;
}

std::string GenChkptVolMol::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "pos=" << pos << ", " <<
      "id=" << id << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "diffusion_time=" << diffusion_time << ", " <<
      "birthday=" << birthday << ", " <<
      "flags=" << flags << ", " <<
      "unimol_rxn_time=" << unimol_rxn_time;
  return ss.str();
}

py::class_<ChkptVolMol> define_pybinding_ChkptVolMol(py::module& m) {
  return py::class_<ChkptVolMol, BaseChkptMol, std::shared_ptr<ChkptVolMol>>(m, "ChkptVolMol", "Class representing a checkpointed volume molecule.\nNot to be used directly.\n")
      .def(
          py::init<
            const Vec3&,
            const int,
            std::shared_ptr<Species>,
            const double,
            const double,
            const int,
            const double
          >(),
          py::arg("pos"),
          py::arg("id"),
          py::arg("species"),
          py::arg("diffusion_time"),
          py::arg("birthday"),
          py::arg("flags"),
          py::arg("unimol_rxn_time") = FLT_UNSET
      )
      .def("check_semantics", &ChkptVolMol::check_semantics)
      .def("__copy__", &ChkptVolMol::copy_chkpt_vol_mol)
      .def("__deepcopy__", &ChkptVolMol::deepcopy_chkpt_vol_mol, py::arg("memo"))
      .def("__str__", &ChkptVolMol::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &ChkptVolMol::__eq__, py::arg("other"))
      .def("dump", &ChkptVolMol::dump)
      .def_property("pos", &ChkptVolMol::get_pos, &ChkptVolMol::set_pos)
    ;
}

std::string GenChkptVolMol::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "chkpt_vol_mol_" + std::to_string(ctx.postinc_counter("chkpt_vol_mol"));
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
  ss << "m.ChkptVolMol(" << nl;
  ss << ind << "id = " << id << "," << nl;
  ss << ind << "species = " << species->export_to_python(out, ctx) << "," << nl;
  ss << ind << "diffusion_time = " << f_to_str(diffusion_time) << "," << nl;
  ss << ind << "birthday = " << f_to_str(birthday) << "," << nl;
  ss << ind << "flags = " << flags << "," << nl;
  if (unimol_rxn_time != FLT_UNSET) {
    ss << ind << "unimol_rxn_time = " << f_to_str(unimol_rxn_time) << "," << nl;
  }
  ss << ind << "pos = " << "m.Vec3(" << f_to_str(pos.x) << ", " << f_to_str(pos.y) << ", " << f_to_str(pos.z)<< ")," << nl;
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

