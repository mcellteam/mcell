
#include "mcell.h"
#include "../generated/gen_release_site.h"
#include "../generated/gen_species.h"

namespace MCell {
namespace API {

void define_binding_Vec3(py::module& m) {
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

  define_binding_Vec3(m);
  define_binding_ReleaseSite(m);
  define_binding_Species(m);

}


} // API
} // MCell


