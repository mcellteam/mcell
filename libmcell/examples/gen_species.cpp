// copyright
// generated info
#include <sstream>

#include "../api/mcell.h"

namespace MCell {
namespace API {


SemRes GenSpecies::check_semantics(std::ostream& out) const {
  if (!is_set(name)) {
    out << get_object_name() << ": Parameter 'name' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenSpecies::to_str() const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "diffusion_constant_2d=" << diffusion_constant_2d  << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d;
  return ss.str();
}


void define_binding_Species(pybind11::module& m) {
    py::class_<MCell::API::Species>(m, "Species")
        .def(
            py::init<
              const std::string&,
              const float_t,
              const float_t
              >(),
            py::arg("name"),
            py::arg("diffusion_constant_2d") = FLT_INVALID,
            py::arg("diffusion_constant_3d") = FLT_INVALID
        )
        .def("check_semantics", &MCell::API::Species::check_semantics_cerr)
        .def("__str__", &MCell::API::Species::to_str)
        .def("dump", &MCell::API::Species::dump)
        .def_property("name", &MCell::API::Species::get_name, &MCell::API::Species::set_name)
        //.def_property("diffusion_constant_2d", &Species::get_site_diameter, &Species::set_site_diameter)
    ;
}

} // namespace API
} // namespace MCell
