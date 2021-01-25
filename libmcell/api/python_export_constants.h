
#ifndef LIBMCELL_API_PYTHON_EXPORT_CONSTANTS_H_
#define LIBMCELL_API_PYTHON_EXPORT_CONSTANTS_H_

#include "datamodel_defines.h"

namespace MCell {
namespace API {

const int FLOAT_OUT_PRECISION = 15; // this is the precision that is used by mdl_to_data_model.py script

const int MAX_SPECIES_NAME_LENGTH = 32;

const char* const PARAMETERS = "parameters";
const char* const SUBSYSTEM = "subsystem";
const char* const GEOMETRY = "geometry";
const char* const INSTANTIATION = "instantiation";
const char* const OBSERVABLES = "observables";
const char* const SIMULATION_STATE = "simulation_state";
const char* const MODEL = "model";

#define CUSTOMIZATION "customization"
const char* const CUSTOM_ARGPARSE_AND_PARAMETERS = "custom_argparse_and_parameters";
const char* const CUSTOM_CONFIG = "custom_config";
const char* const CUSTOM_INIT_AND_RUN = "custom_init_and_run";

const char* const PY_EXT = ".py";
const char* const BNGL_EXT = ".bngl";

const char* const IMPORT = "import";

const char* const MDOT = "m.";

const char* const IND = "    ";
const char* const IND4 = "    "; // this indent must be 4 chars because it is used as indentation in Python
const char* const IND8 = "        ";
const char* const BLOCK_BEGIN1 = "# ---- ";
const char* const BLOCK_BEGIN2 = " ----\n";
const char* const BLOCK_END1 = "# ^^^^ ";
const char* const BLOCK_END2 = " ^^^^\n\n";
const char* const CTOR_END = ")\n\n";

#define PARAM_SEED "SEED"

const char* const PARAM_ITERATIONS = "ITERATIONS";
const char* const PARAM_TIME_STEP = "TIME_STEP";
const char* const PARAM_DUMP = "DUMP";
const char* const PARAM_EXPORT_DATA_MODEL = "EXPORT_DATA_MODEL";

const char* const SET_BNGL_MOLECULE_TYPES_INFO = "set_bngl_molecule_types_info";

const char* const VAR_BNGL_PARAMS = "bngl_params";

const char* const VEC3 = "Vec3";

const char* const COUNT_PREFIX = "count_";
const char* const COUNT_TERM_PREFIX = "cterm_";
const char* const UNNAMED_REACTION_RULE_PREFIX = "unnamed_reaction_rule_";
const char* const SURFACE_CLASS_PREFIX = "surface_class_";
const char* const MOLECULE_LIST_PREFIX = "molecule_list_";


const char* const VIZ_OUTPUT_NAME = "viz_output";
const char* const DEFAULT_VIZ_OUTPUT_FILENAME_PREFIX = "./" VALUE_VIZ_DATA "/seed_' + str(" PARAM_SEED ").zfill(5) + '/Scene";
const char* const DEFAULT_RXN_OUTPUT_FILENAME_PREFIX = "./" VALUE_REACT_DATA "/seed_' + str(" PARAM_SEED ").zfill(5) + '/";

const char* const INTERPRETER = "#!/usr/bin/env python3\n\n";

const char* const GENERATED_WARNING =
  "# WARNING: This is an automatically generated file and will be overwritten\n"
  "#          by CellBlender on the next model export.\n";

const char* const BASE_MODEL_IMPORTS =
    "import sys\n"
    "import os\n"
;

const char* const IMPORT_OS =
    "import os\n"
;

#define MODEL_PATH "MODEL_PATH"

const char* const MODEL_PATH_SETUP =
    MODEL_PATH " = os.path.dirname(os.path.abspath(__file__))\n"
;

const char* const MCELL_PATH_SETUP =
    "# ---- import mcell module located in directory ----\n"
    "# ---- specified by system variable MCELL_PATH  ----\n"
    "MCELL_PATH = os.environ.get('MCELL_PATH', '')\n"
    "if MCELL_PATH:\n"
    "    lib_path = os.path.join(MCELL_PATH, 'lib')\n"
    "    if os.path.exists(os.path.join(lib_path, 'mcell.so')) or \\\n"
    "        os.path.exists(os.path.join(lib_path, 'mcell.pyd')):\n"
    "        sys.path.append(lib_path)\n"
    "    else:\n"
    "        print(\"Error: Python module mcell.so or mcell.pyd was not found in \"\n"
    "              \"directory '\" + lib_path + \"' constructed from system variable \"\n"
    "              \"MCELL_PATH.\")\n"
    "        sys.exit(1)\n"
    "else:\n"
    "    print(\"Error: system variable MCELL_PATH that is used to find the mcell \"\n"
    "          \"library was not set.\")\n"
    "    sys.exit(1)\n"
;

const char* const MCELL_IMPORT = "import mcell as m\n\n";

static std::string get_customization_import(const std::string& customization_module) {
  return
      std::string("if os.path.exists(os.path.join(") + MODEL_PATH + ", '" + customization_module + ".py')):\n"
      "    import " + customization_module + "\n"
      "else:\n"
      "    " + customization_module + " = None\n"
  ;
}

static std::string get_argparse_w_customization_begin(
    const std::string& parameters_module, const std::string& customization_module) {

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
      "        " + parameters_module + ".SEED = int(sys.argv[2])\n"
  ;
}

const char* const CHECKPOINT_ITERATION = "checkpoint_iteration";

static std::string get_argparse_checkpoint_iteration(
    const std::string& parameters_module) {

  return
      std::string("    elif len(sys.argv) == 3 and sys.argv[1] == '-chkpt':\n") +
      "        " + CHECKPOINT_ITERATION + " = int(sys.argv[2])\n" +
      "    elif len(sys.argv) == 5 and sys.argv[1] == '-seed' and sys.argv[3] == '-chkpt':\n" +
      "        " + parameters_module + ".SEED = int(sys.argv[2])\n" +
      "        " + CHECKPOINT_ITERATION + " = int(sys.argv[4])\n"
  ;
}

static std::string get_resume_from_checkpoint_code(
    const std::string& parameters_module) {

  return
      "# resume simulation if a checkpoint was created\n"
      "checkpoint_dir = m.run_utils.get_last_checkpoint_dir(" + parameters_module + ".SEED)\n"
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

static std::string get_argparse_w_customization_end() {
  return
      "    else:\n"
      "        print(\"Error: invalid command line arguments\")\n"
      "        print(\"  usage: \" + sys.argv[0] + \"[-seed N]\")\n"
      "        sys.exit(1)\n";
}


static std::string get_user_defined_configuration(const std::string& customization_module) {
  return
      "if " + customization_module + " and '" + CUSTOM_CONFIG + "' in dir(" + customization_module + "):\n"
      "    # user-defined model configuration\n"
      "    " + customization_module + "." + CUSTOM_CONFIG + "(" + MODEL + ")\n"
  ;
}

static std::string get_abs_path(const std::string file) {
  return std::string("os.path.join(") + MODEL_PATH + ", '" + file + "')";
}

static std::string get_import(const std::string module) {
  return "import " + module + "\n";
}

static std::string get_import_star(const std::string module) {
  return "from " + module + " import *\n";
}

const char* const REGION_ALL_NAME = "ALL";
const char* const REGION_ALL_SUFFIX = "[ALL]";
const char* const NULL_PRODUCTS = "NULL";

const char* const COUNT = "COUNT";
const char* const WORLD = "WORLD";
const char* const WORLD_FIRST_UPPER = "World";

const char* const REV_RXN_SUFFIX = "_rev";

const char* const TEMPLATE_CUSTOM_ARGPARSE_AND_PARAMETERS =
    "\"\"\"\n"
    "def custom_argparse_and_parameters():\n"
    "    # When uncommented, this function is called to parse \n"
    "    # custom commandline arguments.\n"
    "    # It is executed before any of the automatically generated \n"
    "    # parameter values are set so one can override the parameter \n"
    "    # values here as well.\n"
    "    pass\n"
    "\"\"\"\n"
;

const char* const TEMPLATE_CUSTOM_CONFIG =
    "\"\"\"\n"
    "def custom_config(model):\n"
    "    # When uncommented, this function is called to set custom\n"
    "    # model configuration.\n"
    "    # It is executed after basic parameter setup is done and \n"
    "    # before any components are added to the model. \n"
    "    pass\n"
    "\"\"\"\n"
;

const char* const TEMPLATE_CUSTOM_INIT_AND_RUN =
    "\"\"\"\n"
    "def custom_init_and_run(model):\n"
    "    # When uncommented, this function is called after all the model\n"
    "    # components defined in CellBlender were added to the model.\n"
    "    # It allows to add additional model components before initialization \n"
    "    # is done and then to customize how simulation is ran.\n"
    "    model.initialize()\n"
    "    model.run_iterations(parameters.ITERATIONS)\n"
    "    model.end_simulation()\n"
    "\"\"\"\n"
;

} // namespace API
} // namespace MCell

#endif // LIBMCELL_API_PYTHON_EXPORT_CONSTANTS_H_
