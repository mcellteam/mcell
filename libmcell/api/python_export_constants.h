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

#ifndef LIBMCELL_API_PYTHON_EXPORT_CONSTANTS_H_
#define LIBMCELL_API_PYTHON_EXPORT_CONSTANTS_H_

#include "datamodel_defines.h"

namespace MCell {
namespace API {

const int FLOAT_OUT_PRECISION = 15; // this is the precision that is used by mdl_to_data_model.py script

const int MAX_SPECIES_NAME_LENGTH = 32;

const char* const SHARED = "shared";
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

const char* const PARAMETER_OVERRIDES = "parameter_overrides";
const char* const NOT_DEFINED = "not_defined";

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

const char* const COMPONENT_PREFIX = "comp_";
const char* const ELEMENTARY_MOLECULE_TYPE_PREFIX = "emt_";
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

const char* const IMPORT_SYS_OS =
    "import sys\n"
    "import os\n"
;

const char* const IMPORT_OS = "import os\n";
const char* const IMPORT_MATH = "import math\n";
const char* const IMPORT_SHARED = "import shared\n";

#define MODEL_PATH "MODEL_PATH"

const char* const MODEL_PATH_SETUP =
    MODEL_PATH " = os.path.dirname(os.path.abspath(__file__))\n"
;

const char* const MCELL_PATH_SETUP =
    "\n# ---- import mcell module located in directory specified by system variable MCELL_PATH  ----\n\n"
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

const char* const IMPORT_MCELL_AS_M = "import mcell as m\n\n";

const char* const CHECKPOINT_ITERATION = "checkpoint_iteration";

const char* const REGION_ALL_NAME = "ALL";
const char* const REGION_ALL_SUFFIX = "[ALL]";
const char* const NULL_PRODUCTS = "NULL";

const char* const COUNT = "COUNT";
const char* const WORLD = "WORLD";
const char* const WORLD_FIRST_UPPER = "World";

const char* const REV_RXN_SUFFIX = "_rev";

std::string get_customization_import(const std::string& customization_module);
std::string get_argparse_w_customization_begin(const std::string& customization_module);
std::string get_argparse_checkpoint_iteration();
std::string get_resume_from_checkpoint_code();
std::string get_argparse_w_customization_end();
std::string get_user_defined_configuration(const std::string& customization_module);
std::string get_abs_path(const std::string file);
std::string get_import(const std::string module);
std::string get_import_star(const std::string module);

const char* const TEMPLATE_CUSTOM_ARGPARSE_AND_PARAMETERS =
    "\"\"\"\n"
    "def custom_argparse_and_parameters():\n"
    "    # When uncommented, this function is called to parse \n"
    "    # custom commandline arguments.\n"
    "    # It is executed before any of the automatically generated \n"
    "    # parameter values are set so one can override the parameter \n"
    "    # values here as well.\n"
    "    # To override parameter values, add or overwrite an item in dictionary\n"
    "    # shared.parameter_overrides, e.g. shared.parameter_overrides['SEED'] = 10\n"
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

std::string get_template_custom_init_and_run(const std::string& parameters_module);

} // namespace API
} // namespace MCell

#endif // LIBMCELL_API_PYTHON_EXPORT_CONSTANTS_H_
