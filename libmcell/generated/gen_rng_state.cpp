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
#include "gen_rng_state.h"
#include "api/rng_state.h"

namespace MCell {
namespace API {

void GenRngState::check_semantics() const {
  if (!is_set(randcnt)) {
    throw ValueError("Parameter 'randcnt' must be set.");
  }
  if (!is_set(aa)) {
    throw ValueError("Parameter 'aa' must be set.");
  }
  if (!is_set(bb)) {
    throw ValueError("Parameter 'bb' must be set.");
  }
  if (!is_set(cc)) {
    throw ValueError("Parameter 'cc' must be set.");
  }
  if (!is_set(randslr)) {
    throw ValueError("Parameter 'randslr' must be set and the value must not be an empty list.");
  }
  if (!is_set(mm)) {
    throw ValueError("Parameter 'mm' must be set and the value must not be an empty list.");
  }
  if (!is_set(rngblocks)) {
    throw ValueError("Parameter 'rngblocks' must be set.");
  }
}

void GenRngState::set_initialized() {
  initialized = true;
}

void GenRngState::set_all_attributes_as_default_or_unset() {
  class_name = "RngState";
  randcnt = 0;
  aa = 0;
  bb = 0;
  cc = 0;
  randslr = std::vector<uint64_t>();
  mm = std::vector<uint64_t>();
  rngblocks = 0;
}

std::shared_ptr<RngState> GenRngState::copy_rng_state() const {
  std::shared_ptr<RngState> res = std::make_shared<RngState>(DefaultCtorArgType());
  res->class_name = class_name;
  res->randcnt = randcnt;
  res->aa = aa;
  res->bb = bb;
  res->cc = cc;
  res->randslr = randslr;
  res->mm = mm;
  res->rngblocks = rngblocks;

  return res;
}

std::shared_ptr<RngState> GenRngState::deepcopy_rng_state(py::dict) const {
  std::shared_ptr<RngState> res = std::make_shared<RngState>(DefaultCtorArgType());
  res->class_name = class_name;
  res->randcnt = randcnt;
  res->aa = aa;
  res->bb = bb;
  res->cc = cc;
  res->randslr = randslr;
  res->mm = mm;
  res->rngblocks = rngblocks;

  return res;
}

bool GenRngState::__eq__(const RngState& other) const {
  return
    randcnt == other.randcnt &&
    aa == other.aa &&
    bb == other.bb &&
    cc == other.cc &&
    randslr == other.randslr &&
    mm == other.mm &&
    rngblocks == other.rngblocks;
}

bool GenRngState::eq_nonarray_attributes(const RngState& other, const bool ignore_name) const {
  return
    randcnt == other.randcnt &&
    aa == other.aa &&
    bb == other.bb &&
    cc == other.cc &&
    true /*randslr*/ &&
    true /*mm*/ &&
    rngblocks == other.rngblocks;
}

std::string GenRngState::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "randcnt=" << randcnt << ", " <<
      "aa=" << aa << ", " <<
      "bb=" << bb << ", " <<
      "cc=" << cc << ", " <<
      "randslr=" << vec_nonptr_to_str(randslr, all_details, ind + "  ") << ", " <<
      "mm=" << vec_nonptr_to_str(mm, all_details, ind + "  ") << ", " <<
      "rngblocks=" << rngblocks;
  return ss.str();
}

py::class_<RngState> define_pybinding_RngState(py::module& m) {
  return py::class_<RngState, std::shared_ptr<RngState>>(m, "RngState", "Internal checkpointing structure holding state of the random number generator.")
      .def(
          py::init<
            const uint64_t,
            const uint64_t,
            const uint64_t,
            const uint64_t,
            const std::vector<uint64_t>,
            const std::vector<uint64_t>,
            const uint64_t
          >(),
          py::arg("randcnt"),
          py::arg("aa"),
          py::arg("bb"),
          py::arg("cc"),
          py::arg("randslr"),
          py::arg("mm"),
          py::arg("rngblocks")
      )
      .def("check_semantics", &RngState::check_semantics)
      .def("__copy__", &RngState::copy_rng_state)
      .def("__deepcopy__", &RngState::deepcopy_rng_state, py::arg("memo"))
      .def("__str__", &RngState::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &RngState::__eq__, py::arg("other"))
      .def("dump", &RngState::dump)
      .def_property("randcnt", &RngState::get_randcnt, &RngState::set_randcnt)
      .def_property("aa", &RngState::get_aa, &RngState::set_aa)
      .def_property("bb", &RngState::get_bb, &RngState::set_bb)
      .def_property("cc", &RngState::get_cc, &RngState::set_cc)
      .def_property("randslr", &RngState::get_randslr, &RngState::set_randslr, py::return_value_policy::reference, "Must contain RNG_SIZE items.")
      .def_property("mm", &RngState::get_mm, &RngState::set_mm, py::return_value_policy::reference, "Must contain RNG_SIZE items.")
      .def_property("rngblocks", &RngState::get_rngblocks, &RngState::set_rngblocks, "Must contain RNG_SIZE items.")
    ;
}

std::string GenRngState::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "rng_state_" + std::to_string(ctx.postinc_counter("rng_state"));
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
  ss << "m.RngState(" << nl;
  ss << ind << "randcnt = " << randcnt << "," << nl;
  ss << ind << "aa = " << aa << "," << nl;
  ss << ind << "bb = " << bb << "," << nl;
  ss << ind << "cc = " << cc << "," << nl;
  ss << ind << "randslr = " << export_vec_randslr(out, ctx, exported_name) << "," << nl;
  ss << ind << "mm = " << export_vec_mm(out, ctx, exported_name) << "," << nl;
  ss << ind << "rngblocks = " << rngblocks << "," << nl;
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

std::string GenRngState::export_vec_randslr(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < randslr.size(); i++) {
    const auto& item = randslr[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << item << ", ";
  }
  ss << "]";
  return ss.str();
}

std::string GenRngState::export_vec_mm(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < mm.size(); i++) {
    const auto& item = mm[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << item << ", ";
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

