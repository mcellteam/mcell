/******************************************************************************
 *
 * Copyright (C) 2021 by
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

#include "api/python_export_constants.h"

namespace MCell {
namespace API {

std::string get_customization_import(const std::string& customization_module) {
  return
      std::string("if os.path.exists(os.path.join(") + MODEL_PATH + ", '" + customization_module + ".py')):\n"
      "    import " + customization_module + "\n"
      "else:\n"
      "    " + customization_module + " = None\n"
  ;
}

std::string get_argparse_w_customization_begin(const std::string& customization_module) {

  return
      "if " + customization_module + " and '" + CUSTOM_ARGPARSE_AND_PARAMETERS + "' in dir(" + customization_module + "):\n"
      "    # custom argument processing and parameter setup\n"
      "    " + customization_module + "." + CUSTOM_ARGPARSE_AND_PARAMETERS + "()\n"
      "else:\n"
      "    if len(sys.argv) == 1:\n"
      "        # no arguments\n"
      "        pass\n"
      "    elif len(sys.argv) == 3 and sys.argv[1] == '-seed':\n"
      "        # overwrite value of seed defined in module parameters\n"
      "        " + SHARED + "." + PARAMETER_OVERRIDES + "['" + PARAM_SEED + "'] = int(sys.argv[2])\n"
  ;
}

std::string get_argparse_checkpoint_iteration() {

  return
      std::string("    elif len(sys.argv) == 3 and sys.argv[1] == '-chkpt':\n") +
      "        " + CHECKPOINT_ITERATION + " = int(sys.argv[2])\n" +
      "    elif len(sys.argv) == 5 and sys.argv[1] == '-seed' and sys.argv[3] == '-chkpt':\n" +
      "        " + SHARED + "." + PARAMETER_OVERRIDES + "['" + PARAM_SEED + "'] = int(sys.argv[2])\n" +
      "        " + CHECKPOINT_ITERATION + " = int(sys.argv[4])\n"
  ;
}

std::string get_resume_from_checkpoint_code() {

  return
      "# resume simulation if a checkpoint was created\n"
      "checkpoint_dir = m.run_utils.get_last_checkpoint_dir(" PARAM_SEED ")\n"
      "if checkpoint_dir:\n"
      "    # change sys.path so that the only the checkpointed files are loaded\n"
      "    sys.path = m.run_utils.remove_cwd(sys.path)\n"
      "    sys.path.append(checkpoint_dir)\n"
      "    \n"
      "    # prepare import of the 'model' module from the checkpoint\n"
      "    model_spec = importlib.util.spec_from_file_location(\n"
      "        'model', os.path.join(checkpoint_dir, 'model.py'))\n"
      "    model_module = importlib.util.module_from_spec(model_spec)\n"
      "    \n"
      "    # run import, this also resumes simulation from the checkpoint\n"
      "    model_spec.loader.exec_module(model_module)\n"
      "\n"
      "    # exit after simulation has finished successfully\n"
      "    sys.exit(0)\n\n"
 ;
}

std::string get_argparse_w_customization_end() {
  return
      "    else:\n"
      "        print(\"Error: invalid command line arguments\")\n"
      "        print(\"  usage: \" + sys.argv[0] + \"[-seed N]\")\n"
      "        sys.exit(1)\n";
}


std::string get_user_defined_configuration(const std::string& customization_module) {
  return
      "if " + customization_module + " and '" + CUSTOM_CONFIG + "' in dir(" + customization_module + "):\n"
      "    # user-defined model configuration\n"
      "    " + customization_module + "." + CUSTOM_CONFIG + "(" + MODEL + ")\n"
  ;
}

std::string get_abs_path(const std::string file) {
  return std::string("os.path.join(") + MODEL_PATH + ", '" + file + "')";
}

std::string get_import(const std::string module) {
  return "import " + module + "\n";
}

std::string get_import_star(const std::string module) {
  return "from " + module + " import *\n";
}


std::string get_template_custom_init_and_run(const std::string& parameters_module) {
  return
      "\"\"\"\n"
      "def custom_init_and_run(model):\n"
      "    # When uncommented, this function is called after all the model\n"
      "    # components defined in CellBlender were added to the model.\n"
      "    # It allows to add additional model components before initialization \n"
      "    # is done and then to customize how simulation is ran.\n"
      "    # The module parameters must be imported locally otherwise "
      "    # changes to shared.parameter_overrides done elsewhere won't be applied.\n"
      "    import " + parameters_module + " as " + PARAMETERS + "\n" +
      "    model.initialize()\n"
      "    model.run_iterations(parameters.ITERATIONS)\n"
      "    model.end_simulation()\n"
      "\"\"\"\n";
}

} // namespace API
} // namespace MCell

