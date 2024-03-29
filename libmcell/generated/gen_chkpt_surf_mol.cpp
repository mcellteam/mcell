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
#include "gen_chkpt_surf_mol.h"
#include "api/chkpt_surf_mol.h"
#include "api/geometry_object.h"
#include "api/species.h"

namespace MCell {
namespace API {

void GenChkptSurfMol::check_semantics() const {
  if (!is_set(pos)) {
    throw ValueError("Parameter 'pos' must be set.");
  }
  if (!is_set(orientation)) {
    throw ValueError("Parameter 'orientation' must be set.");
  }
  if (!is_set(geometry_object)) {
    throw ValueError("Parameter 'geometry_object' must be set.");
  }
  if (!is_set(wall_index)) {
    throw ValueError("Parameter 'wall_index' must be set.");
  }
  if (!is_set(grid_tile_index)) {
    throw ValueError("Parameter 'grid_tile_index' must be set.");
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

void GenChkptSurfMol::set_initialized() {
  if (is_set(geometry_object)) {
    geometry_object->set_initialized();
  }
  if (is_set(species)) {
    species->set_initialized();
  }
  initialized = true;
}

void GenChkptSurfMol::set_all_attributes_as_default_or_unset() {
  class_name = "ChkptSurfMol";
  pos = VEC2_UNSET;
  orientation = Orientation::NOT_SET;
  geometry_object = nullptr;
  wall_index = INT_UNSET;
  grid_tile_index = INT_UNSET;
  id = INT_UNSET;
  species = nullptr;
  diffusion_time = FLT_UNSET;
  birthday = FLT_UNSET;
  flags = INT_UNSET;
  unimol_rxn_time = FLT_UNSET;
}

std::shared_ptr<ChkptSurfMol> GenChkptSurfMol::copy_chkpt_surf_mol() const {
  std::shared_ptr<ChkptSurfMol> res = std::make_shared<ChkptSurfMol>(DefaultCtorArgType());
  res->class_name = class_name;
  res->pos = pos;
  res->orientation = orientation;
  res->geometry_object = geometry_object;
  res->wall_index = wall_index;
  res->grid_tile_index = grid_tile_index;
  res->id = id;
  res->species = species;
  res->diffusion_time = diffusion_time;
  res->birthday = birthday;
  res->flags = flags;
  res->unimol_rxn_time = unimol_rxn_time;

  return res;
}

std::shared_ptr<ChkptSurfMol> GenChkptSurfMol::deepcopy_chkpt_surf_mol(py::dict) const {
  std::shared_ptr<ChkptSurfMol> res = std::make_shared<ChkptSurfMol>(DefaultCtorArgType());
  res->class_name = class_name;
  res->pos = pos;
  res->orientation = orientation;
  res->geometry_object = is_set(geometry_object) ? geometry_object->deepcopy_geometry_object() : nullptr;
  res->wall_index = wall_index;
  res->grid_tile_index = grid_tile_index;
  res->id = id;
  res->species = is_set(species) ? species->deepcopy_species() : nullptr;
  res->diffusion_time = diffusion_time;
  res->birthday = birthday;
  res->flags = flags;
  res->unimol_rxn_time = unimol_rxn_time;

  return res;
}

bool GenChkptSurfMol::__eq__(const ChkptSurfMol& other) const {
  return
    pos == other.pos &&
    orientation == other.orientation &&
    (
      (is_set(geometry_object)) ?
        (is_set(other.geometry_object) ?
          (geometry_object->__eq__(*other.geometry_object)) : 
          false
        ) :
        (is_set(other.geometry_object) ?
          false :
          true
        )
     )  &&
    wall_index == other.wall_index &&
    grid_tile_index == other.grid_tile_index &&
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

bool GenChkptSurfMol::eq_nonarray_attributes(const ChkptSurfMol& other, const bool ignore_name) const {
  return
    pos == other.pos &&
    orientation == other.orientation &&
    (
      (is_set(geometry_object)) ?
        (is_set(other.geometry_object) ?
          (geometry_object->__eq__(*other.geometry_object)) : 
          false
        ) :
        (is_set(other.geometry_object) ?
          false :
          true
        )
     )  &&
    wall_index == other.wall_index &&
    grid_tile_index == other.grid_tile_index &&
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

std::string GenChkptSurfMol::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "pos=" << pos << ", " <<
      "orientation=" << orientation << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "grid_tile_index=" << grid_tile_index << ", " <<
      "id=" << id << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "diffusion_time=" << diffusion_time << ", " <<
      "birthday=" << birthday << ", " <<
      "flags=" << flags << ", " <<
      "unimol_rxn_time=" << unimol_rxn_time;
  return ss.str();
}

py::class_<ChkptSurfMol> define_pybinding_ChkptSurfMol(py::module& m) {
  return py::class_<ChkptSurfMol, BaseChkptMol, std::shared_ptr<ChkptSurfMol>>(m, "ChkptSurfMol", "Class representing a checkpointed surface molecule.\nNot to be used directly.\n")
      .def(
          py::init<
            const Vec2&,
            const Orientation,
            std::shared_ptr<GeometryObject>,
            const int,
            const int,
            const int,
            std::shared_ptr<Species>,
            const double,
            const double,
            const int,
            const double
          >(),
          py::arg("pos"),
          py::arg("orientation"),
          py::arg("geometry_object"),
          py::arg("wall_index"),
          py::arg("grid_tile_index"),
          py::arg("id"),
          py::arg("species"),
          py::arg("diffusion_time"),
          py::arg("birthday"),
          py::arg("flags"),
          py::arg("unimol_rxn_time") = FLT_UNSET
      )
      .def("check_semantics", &ChkptSurfMol::check_semantics)
      .def("__copy__", &ChkptSurfMol::copy_chkpt_surf_mol)
      .def("__deepcopy__", &ChkptSurfMol::deepcopy_chkpt_surf_mol, py::arg("memo"))
      .def("__str__", &ChkptSurfMol::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &ChkptSurfMol::__eq__, py::arg("other"))
      .def("dump", &ChkptSurfMol::dump)
      .def_property("pos", &ChkptSurfMol::get_pos, &ChkptSurfMol::set_pos)
      .def_property("orientation", &ChkptSurfMol::get_orientation, &ChkptSurfMol::set_orientation)
      .def_property("geometry_object", &ChkptSurfMol::get_geometry_object, &ChkptSurfMol::set_geometry_object)
      .def_property("wall_index", &ChkptSurfMol::get_wall_index, &ChkptSurfMol::set_wall_index)
      .def_property("grid_tile_index", &ChkptSurfMol::get_grid_tile_index, &ChkptSurfMol::set_grid_tile_index)
    ;
}

std::string GenChkptSurfMol::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "chkpt_surf_mol_" + std::to_string(ctx.postinc_counter("chkpt_surf_mol"));
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
  ss << "m.ChkptSurfMol(" << nl;
  ss << ind << "id = " << id << "," << nl;
  ss << ind << "species = " << species->export_to_python(out, ctx) << "," << nl;
  ss << ind << "diffusion_time = " << f_to_str(diffusion_time) << "," << nl;
  ss << ind << "birthday = " << f_to_str(birthday) << "," << nl;
  ss << ind << "flags = " << flags << "," << nl;
  if (unimol_rxn_time != FLT_UNSET) {
    ss << ind << "unimol_rxn_time = " << f_to_str(unimol_rxn_time) << "," << nl;
  }
  ss << ind << "pos = " << "m.Vec2(" << f_to_str(pos.u) << ", " << f_to_str(pos.v)<< ")," << nl;
  ss << ind << "orientation = " << orientation << "," << nl;
  ss << ind << "geometry_object = " << geometry_object->export_to_python(out, ctx) << "," << nl;
  ss << ind << "wall_index = " << wall_index << "," << nl;
  ss << ind << "grid_tile_index = " << grid_tile_index << "," << nl;
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

