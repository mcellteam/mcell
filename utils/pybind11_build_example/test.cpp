
#define _hypot hypot
#include <cmath>
#include <pybind11/pybind11.h>

#include <string>
#include <vector>

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

enum class Orientation {
  DOWN = -1,
  NONE = 0,
  UP = 1,
  NOT_SET = 2,
  ANY = 3,
  DEFAULT = 4
};

#define POS_INVALID 0

struct Vec3 {
  Vec3() = default;
  Vec3(const Vec3& a) { }
  Vec3(const float x_, const float y_, const float z_) { x = x_; y = y_; z = z_; }
  Vec3(const float xyz) { x = xyz; y = xyz; z = xyz; }
  Vec3(const std::vector<float>& xyz) { assert(xyz.size() == 3); x = xyz[0]; y = xyz[1]; z = xyz[2]; }

  bool is_valid() const { return !(x == POS_INVALID || y == POS_INVALID || z == POS_INVALID); }

  std::string to_string() const;
  void dump(const std::string extra_comment, const std::string ind) const;

  float x, y, z;
};


void define_pybinding_Vec3(py::module& m) {
  py::class_<Vec3>(m, "Vec3")
      .def(
          py::init<>()
      )
      .def(
          py::init<const float>(),
          py::arg("xyz")
      )
      .def(
          py::init<const float, const float, const float>(),
          py::arg("x"), py::arg("y"), py::arg("z")
      )
      .def("__add__", [](const Vec3& a, const Vec3& b) { return Vec3(a.x + b.x); } )
      .def("__sub__", [](const Vec3& a, const Vec3& b) { return Vec3(a.x - b.x); } )
      .def("__mul__", [](const Vec3& a, const Vec3& b) { return Vec3(a.x * b.x); } )
      .def("__truediv__", [](const Vec3& a, const Vec3& b) { return Vec3(a.x / b.x); } )
      .def("__eq__",  [](const Vec3& a, const Vec3& b) { return a.x == b.x; } )
      .def("__str__",  [](const Vec3& a)
          { return "(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + std::to_string(a.z) + ")"; } )
      .def("__repr__",  [](const Vec3& a)
          { return "(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + std::to_string(a.z) + ")"; } )
      .def("to_list",  [](const Vec3& a) { return std::vector<float>{a.x, a.y, a.z}; } )
      .def_readwrite("x", &Vec3::x)
      .def_readwrite("y", &Vec3::y)
      .def_readwrite("z", &Vec3::z)
  ;
}


PYBIND11_MODULE(mcell, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    py::enum_<Orientation>(m, "Orientation", py::arithmetic())
      .value("DOWN", Orientation::DOWN)
      .value("NONE", Orientation::NONE)
      .value("UP", Orientation::UP)
      .value("NOT_SET", Orientation::NOT_SET)
      .value("ANY", Orientation::ANY)
      .value("DEFAULT", Orientation::DEFAULT)
      .export_values();

    define_pybinding_Vec3(m);

    m.def("add", &add, "A function which adds two numbers");
}
