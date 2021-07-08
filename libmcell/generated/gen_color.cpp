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
#include "gen_color.h"
#include "api/color.h"

namespace MCell {
namespace API {

void GenColor::check_semantics() const {
}

void GenColor::set_initialized() {
  initialized = true;
}

void GenColor::set_all_attributes_as_default_or_unset() {
  class_name = "Color";
  red = FLT_UNSET;
  green = FLT_UNSET;
  blue = FLT_UNSET;
  alpha = 1;
  rgba = 0;
}

std::shared_ptr<Color> GenColor::copy_color() const {
  std::shared_ptr<Color> res = std::make_shared<Color>(DefaultCtorArgType());
  res->class_name = class_name;
  res->red = red;
  res->green = green;
  res->blue = blue;
  res->alpha = alpha;
  res->rgba = rgba;

  return res;
}

std::shared_ptr<Color> GenColor::deepcopy_color(py::dict) const {
  std::shared_ptr<Color> res = std::make_shared<Color>(DefaultCtorArgType());
  res->class_name = class_name;
  res->red = red;
  res->green = green;
  res->blue = blue;
  res->alpha = alpha;
  res->rgba = rgba;

  return res;
}

bool GenColor::__eq__(const Color& other) const {
  return
    red == other.red &&
    green == other.green &&
    blue == other.blue &&
    alpha == other.alpha &&
    rgba == other.rgba;
}

bool GenColor::eq_nonarray_attributes(const Color& other, const bool ignore_name) const {
  return
    red == other.red &&
    green == other.green &&
    blue == other.blue &&
    alpha == other.alpha &&
    rgba == other.rgba;
}

std::string GenColor::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "red=" << red << ", " <<
      "green=" << green << ", " <<
      "blue=" << blue << ", " <<
      "alpha=" << alpha << ", " <<
      "rgba=" << rgba;
  return ss.str();
}

py::class_<Color> define_pybinding_Color(py::module& m) {
  return py::class_<Color, std::shared_ptr<Color>>(m, "Color", "Represents color with alpha component.\nProvides two means to set value, either red, green, blue and alpha, \nor rgba. If both color individual components and rgba are set in initialization,\nthe individual components are used.\n \n")
      .def(
          py::init<
            const double,
            const double,
            const double,
            const double,
            const uint
          >(),
          py::arg("red") = FLT_UNSET,
          py::arg("green") = FLT_UNSET,
          py::arg("blue") = FLT_UNSET,
          py::arg("alpha") = 1,
          py::arg("rgba") = 0
      )
      .def("check_semantics", &Color::check_semantics)
      .def("__copy__", &Color::copy_color)
      .def("__deepcopy__", &Color::deepcopy_color, py::arg("memo"))
      .def("__str__", &Color::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Color::__eq__, py::arg("other"))
      .def("dump", &Color::dump)
      .def_property("red", &Color::get_red, &Color::set_red, "Red component in range 0-1.")
      .def_property("green", &Color::get_green, &Color::set_green, "Green component in range 0-1.")
      .def_property("blue", &Color::get_blue, &Color::set_blue, "Blue component in range 0-1.")
      .def_property("alpha", &Color::get_alpha, &Color::set_alpha, "Alpha component in range 0-1. 1 means nontransparent.")
      .def_property("rgba", &Color::get_rgba, &Color::set_rgba, "This attribute provides an alternative way of defining colors by supplying a \n32-bit unsigned integer representation of the color with an aplha channel. \nIn hexadecimal notation the first 2 digits are value for red, second 2 digits are \ngreen, third 2 digits are blue and the last two digits are alpha. \nThe range for each component is thus 0x0-0xFF (0-255). \nExample: 0x0000ffcc represents the same color as rgba(0, 0, 100%, 80%).\nAll values are valid.\n   \n")
    ;
}

std::string GenColor::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "color_" + std::to_string(ctx.postinc_counter("color"));
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
  ss << "m.Color(" << nl;
  if (red != FLT_UNSET) {
    ss << ind << "red = " << f_to_str(red) << "," << nl;
  }
  if (green != FLT_UNSET) {
    ss << ind << "green = " << f_to_str(green) << "," << nl;
  }
  if (blue != FLT_UNSET) {
    ss << ind << "blue = " << f_to_str(blue) << "," << nl;
  }
  if (alpha != 1) {
    ss << ind << "alpha = " << f_to_str(alpha) << "," << nl;
  }
  if (rgba != 0) {
    ss << ind << "rgba = " << rgba << "," << nl;
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

