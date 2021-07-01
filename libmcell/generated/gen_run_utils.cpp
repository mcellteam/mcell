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
#include "gen_run_utils.h"

namespace MCell {
namespace API {

void define_pybinding_run_utils(py::module& m) {
  m.def_submodule("run_utils")
      .def("get_last_checkpoint_dir", &run_utils::get_last_checkpoint_dir, py::arg("seed"), "Searches the directory checkpoints for the last checkpoint for the given \nparameters and returns the directory name if such a directory exists. \nReturns empty string if no checkpoint directory was found.\nCurrently supports only the seed argument.\n\n- seed\n")
      .def("remove_cwd", &run_utils::remove_cwd, py::arg("paths"), "Removes all directory names items pointing to the current working directory from a list and \nreturns a new list.\n\n- paths\n")
    ;
}

} // namespace API
} // namespace MCell

