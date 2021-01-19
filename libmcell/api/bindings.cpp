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

#include "mcell.h"

#include "generated/gen_component_type.h"
#include "generated/gen_component.h"
#include "generated/gen_elementary_molecule_type.h"
#include "generated/gen_elementary_molecule.h"
#include "generated/gen_species.h"
#include "generated/gen_surface_class.h"
#include "generated/gen_reaction_rule.h"
#include "generated/gen_subsystem.h"

#include "generated/gen_region.h"
#include "generated/gen_initial_surface_release.cpp"
#include "generated/gen_surface_region.h"
#include "generated/gen_geometry_object.h"
#include "generated/gen_release_pattern.h"
#include "generated/gen_molecule_release_info.h"
#include "generated/gen_release_site.h"
#include "generated/gen_instantiation.h"

#include "generated/gen_count_term.h"
#include "generated/gen_count.h"
#include "generated/gen_viz_output.h"
#include "generated/gen_observables.h"

#include "generated/gen_config.h"
#include "generated/gen_notifications.h"
#include "generated/gen_warnings.h"
#include "generated/gen_model.h"

#include "generated/gen_region.h"

#include "generated/gen_molecule.h"
#include "generated/gen_wall.h"
#include "generated/gen_introspection.h"

#include "generated/gen_mol_wall_hit_info.h"

#include "generated/gen_wall_wall_hit_info.h"
#include "generated/gen_reaction_info.h"

#include "generated/gen_chkpt_vol_mol.h"
#include "generated/gen_chkpt_surf_mol.h"
#include "generated/gen_rng_state.h"

#include "generated/gen_bngl_utils.h"
#include "generated/gen_geometry_utils.h"
#include "generated/gen_run_utils.h"


#include "generated/gen_constants.h"

#include "world.h"
#include "version.h"


#if __cplusplus < 201402L && !defined(_MSC_VER)
#error "Pybind11 overload requires at least C++14"
#endif

#ifndef PYBIND11_OVERLOAD_CAST
#error "PYBIND11_OVERLOAD_CAST must be defined"
#endif

namespace MCell {
namespace API {


void define_pybinding_Vec3(py::module& m) {
  py::class_<MCell::Vec3>(m, "Vec3")
      .def(
          py::init<>()
      )
      .def(
          py::init<const float_t>(),
          py::arg("xyz")
      )
      .def(
          py::init<const float_t, const float_t, const float_t>(),
          py::arg("x"), py::arg("y"), py::arg("z")
      )
      .def("__add__", [](const Vec3& a, const Vec3& b) { return Vec3(a + b); } )
      .def("__sub__", [](const Vec3& a, const Vec3& b) { return Vec3(a - b); } )
      .def("__mul__", [](const Vec3& a, const Vec3& b) { return Vec3(a * b); } )
      .def("__truediv__", [](const Vec3& a, const Vec3& b) { return Vec3(a / b); } )
      .def("__eq__",  [](const Vec3& a, const Vec3& b) { return a == b; } )
      .def("__str__",  [](const Vec3& a)
          { return "(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + std::to_string(a.z) + ")"; } )
      .def("__repr__",  [](const Vec3& a)
          { return "(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + std::to_string(a.z) + ")"; } )
      .def("to_list",  [](const Vec3& a) { return std::vector<float_t>{a.x, a.y, a.z}; } )
      .def_readwrite("x", &Vec3::x)
      .def_readwrite("y", &Vec3::y)
      .def_readwrite("z", &Vec3::z)
  ;
}

void define_pybinding_Vec2(py::module& m) {
  py::class_<MCell::Vec2>(m, "Vec2")
      .def(
          py::init<>()
      )
      .def(
          py::init<const float_t>(),
          py::arg("xy")
      )
      .def(
          py::init<const float_t, const float_t>(),
          py::arg("x"), py::arg("y")
      )
      .def("__add__", [](const Vec2& a, const Vec2& b) { return Vec2(a + b); } )
      .def("__sub__", [](const Vec2& a, const Vec2& b) { return Vec2(a - b); } )
      .def("__mul__", [](const Vec2& a, const Vec2& b) { return Vec2(a * b); } )
      .def("__truediv__", [](const Vec2& a, const Vec2& b) { return Vec2(a / b); } )
      .def("__eq__",  [](const Vec2& a, const Vec2& b) { return a == b; } )
      .def("__str__",  [](const Vec2& a)
          { return "(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ")"; } )
      .def("__repr__",  [](const Vec2& a)
          { return "(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ")"; } )
      .def("to_list",  [](const Vec2& a) { return std::vector<float_t>{a.x, a.y}; } )
      .def_readwrite("x", &Vec2::x)
      .def_readwrite("y", &Vec2::y)
      .def_readwrite("u", &Vec2::u)
      .def_readwrite("v", &Vec2::v)
  ;
}

void define_pybinding_IVec3(py::module& m) {
  py::class_<MCell::IVec3>(m, "IVec3")
      .def(
          py::init<>()
      )
      .def(
          py::init<const float_t>(),
          py::arg("xyz")
      )
      .def(
          py::init<const int, const int, const int>(),
          py::arg("x"), py::arg("y"), py::arg("z")
      )
      // FIXME: enabling these operator breaks debug build for unknown reason,
      // release is build ok, one of tests to reproduce:
      // mdl_data_model_pymcell4/tests/mdl/0051_largest_partition_size
      // reenable parts of testing in mcell_tests/tests/pymcell4_positive/0050_vec_operators
      // once fixed
      /*.def(
          py::init<const float_t, const float_t, const float_t>(),
          py::arg("x"), py::arg("y"), py::arg("z")
      )
      .def("__add__", [](const IVec3& a, const IVec3& b) { return IVec3(a + b); } )
      .def("__sub__", [](const IVec3& a, const IVec3& b) { return IVec3(a - b); } )
      .def("__mul__", [](const IVec3& a, const IVec3& b) { return IVec3(a * b); } )
      .def("__truediv__", [](const IVec3& a, const IVec3& b) { return IVec3(a / b); } )
      .def("__eq__",  [](const IVec3& a, const IVec3& b) { return a == b; } )
      */
      .def_readwrite("x", &IVec3::x)
      .def_readwrite("y", &IVec3::y)
      .def_readwrite("z", &IVec3::z)
  ;
}

// all define_binding_* functions are called here
PYBIND11_MODULE(mcell, m) {

  // version
  m.attr("__version__") = py::str(MCELL_VERSION);

  // some classes use enums, must be defined first
  define_pybinding_enums(m);

  define_pybinding_Vec3(m);
  define_pybinding_Vec2(m);
  define_pybinding_IVec3(m);

  define_pybinding_ComponentType(m);
  define_pybinding_Component(m);
  define_pybinding_ElementaryMoleculeType(m);
  define_pybinding_ElementaryMolecule(m);
  define_pybinding_Complex(m);
  define_pybinding_Species(m);
  define_pybinding_SurfaceProperty(m);
  define_pybinding_SurfaceClass(m);
  define_pybinding_ReactionRule(m);
  define_pybinding_Subsystem(m);

  define_pybinding_ReleasePattern(m);
  define_pybinding_MoleculeReleaseInfo(m);
  define_pybinding_ReleaseSite(m);

  define_pybinding_Region(m);
  define_pybinding_InitialSurfaceRelease(m);
  define_pybinding_SurfaceRegion(m);
  define_pybinding_GeometryObject(m);
  define_pybinding_Instantiation(m);

  define_pybinding_CountTerm(m);
  define_pybinding_Count(m);
  define_pybinding_VizOutput(m);
  define_pybinding_Observables(m);

  define_pybinding_Config(m);
  define_pybinding_Notifications(m);
  define_pybinding_Warnings(m);
  define_pybinding_Model(m);

  define_pybinding_Molecule(m);
  define_pybinding_Wall(m);
  define_pybinding_WallWallHitInfo(m);

  define_pybinding_MolWallHitInfo(m);
  define_pybinding_ReactionInfo(m);

  define_pybinding_BaseChkptMol(m);
  define_pybinding_ChkptVolMol(m);
  define_pybinding_ChkptSurfMol(m);
  define_pybinding_RngState(m);

  // constants may reference existing types, must be "bound" later
  define_pybinding_constants(m);

  define_pybinding_bngl_utils(m);
  define_pybinding_geometry_utils(m);
  define_pybinding_run_utils(m);
}


void check_ctrl_c(const float_t current_time, World* world) {
  // make sure to re-acquire lock, PyErr_CheckSignals segfaults otherwise
  py::gil_scoped_acquire acquire;
  if (PyErr_CheckSignals() != 0) {
    std::cout << "Caught Ctrl-C signal in iteration " << current_time << ".\n";
    std::cout << "Flushing buffers.\n";
    release_assert(world != nullptr);
    world->flush_buffers();

    throw py::error_already_set();
  }
}


} // API
} // MCell


