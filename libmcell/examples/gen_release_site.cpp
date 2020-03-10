
#include <sstream>
#include "../api/mcell.h"

namespace MCell {
namespace API {


SemRes GenReleaseSite::check_semantics(std::ostream& out) const {
  if (!is_set(name)) {
    out << get_object_name() << ": Parameter 'name' must be set.\n";
    return SemRes::ERROR;
  }
  if (!is_set(shape)) {
    out << get_object_name() << ": Parameter 'shape' must be set.\n";
    return SemRes::ERROR;
  }
  if (!is_set(molecule)) {
    out << get_object_name() << ": Parameter 'molecule' must be set.\n";
    return SemRes::ERROR;
  }
  return SemRes::OK;
}

std::string GenReleaseSite::to_str() const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "shape = " << shape  << ", " <<
      "location = " << location  << ", " <<
      "molecule = (" << ((molecule != nullptr) ? molecule->to_str() : "null" ) << ")"
      ;
  return ss.str();
}

void define_binding_ReleaseSite(py::module& m) {
    py::class_<ReleaseSite>(m, "ReleaseSite")
        .def(
            py::init<
              const std::string&,
              const std::string&,
              const Species*,
              const Vec3&,
              const float_t,
              const float_t,
              const float_t
              >(),
            py::arg("name"),
            py::arg("shape"),
            py::arg("molecule"),
            py::arg("location") = FLT_INVALID,
            py::arg("site_diameter") = FLT_INVALID,
            py::arg("site_radius") = FLT_INVALID,
            py::arg("release_probability") = FLT_INVALID
        )
        .def("check_semantics", &ReleaseSite::check_semantics_cerr)
        .def("__str__", &MCell::API::ReleaseSite::to_str)
        .def("dump", &MCell::API::ReleaseSite::dump)
        .def_property("name", &ReleaseSite::get_name, &ReleaseSite::set_name)
        .def_property("site_diameter", &ReleaseSite::get_site_diameter, &ReleaseSite::set_site_diameter)
    ;
}

} // namespace API
} // namespace MCell
