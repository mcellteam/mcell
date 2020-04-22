/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#include "../generated/gen_component_type.h"
#include "../generated/gen_component_instance.h"
#include "../generated/gen_molecule_type.h"
#include "../generated/gen_molecule_instance.h"

#include "../generated/gen_species.h"
#include "../generated/gen_release_site.h"
#include "../generated/gen_geometry_object.h"

#include "../generated/gen_subsystem.h"
#include "../generated/gen_instantiation_data.h"
#include "../generated/gen_model.h"

#include "../generated/gen_constants.h"

#if __cplusplus < 201402L
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
  ;
}

// all define_binding_* functions are called here
PYBIND11_MODULE(mcell, m) {

  py::enum_<SemRes>(m, "SemRes", py::arithmetic())
      .value("UNCHECKED", SemRes::UNCHECKED)
      .value("OK", SemRes::OK)
      .value("MESSAGE", SemRes::MESSAGE)
      .value("WARNING", SemRes::WARNING)
      .value("ERROR", SemRes::ERROR)
      .value("FATAL_ERROR", SemRes::FATAL_ERROR)
      .export_values();

  define_pybinding_constants(m);

  define_pybinding_Vec3(m);

  define_pybinding_ComponentType(m);
  define_pybinding_ComponentInstance(m);
  define_pybinding_MoleculeType(m);
  define_pybinding_MoleculeInstance(m);
  define_pybinding_ComplexInstance(m);
  define_pybinding_Species(m);

  define_pybinding_ReleaseSite(m);
  define_pybinding_GeometryObject(m);

  define_pybinding_Subsystem(m);
  define_pybinding_InstantiationData(m);
  define_pybinding_Model(m);

}


} // API
} // MCell


