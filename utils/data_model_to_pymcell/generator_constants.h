
#ifndef SRC4_PYMCELLCONVERTER_GENERATOR_CONSTANTSS_H_
#define SRC4_PYMCELLCONVERTER_GENERATOR_CONSTANTSS_H_

#include "datamodel_defines.h"

namespace MCell {

const int FLOAT_OUT_PRECISION = 15; // this is the precision that is used by mdl_to_data_model.py script

const char* const PARAMETERS = "parameters";
const char* const SUBSYSTEM = "subsystem";
const char* const GEOMETRY = "geometry";
const char* const INSTANTIATION = "instantiation";
const char* const OBSERVABLES = "observables";
const char* const MODEL = "model";

const char* const PY_EXT = ".py";

const char* const IMPORT = "import";

const char* const MDOT = "m.";

const char* const IND = "    ";
const char* const BLOCK_BEGIN1 = "# ---- ";
const char* const BLOCK_BEGIN2 = " ----\n";
const char* const BLOCK_END1 = "# ^^^^ ";
const char* const BLOCK_END2 = " ^^^^\n\n";
const char* const CTOR_END = ")\n\n";

#define PARAM_SEED "SEED_MODULE_LOCAL"
#define FUNCTION_GET_SEED "get_seed"
#define FUNCTION_UPDATE_SEED "update_seed"
const char* const PARAM_ITERATIONS = "ITERATIONS";
const char* const PARAM_TIME_STEP = "TIME_STEP";
const char* const PARAM_DUMP = "DUMP";
const char* const PARAM_EXPORT_DATA_MODEL = "EXPORT_DATA_MODEL";

const char* const VEC3 = "Vec3";

const char* const COUNT_PREFIX = "count_";
const char* const COUNT_TERM_PREFIX = "cterm_";
const char* const UNNAMED_REACTION_RULE_PREFIX = "unnamed_reaction_rule_";
const char* const SURFACE_CLASS_PREFIX = "surface_class_";
const char* const MOLECULE_LIST_PREFIX = "molecule_list_";


const char* const VIZ_OUTPUT_NAME = "viz_output";
const char* const DEFAULT_VIZ_OUTPUT_FILENAME_PREFIX = "./" VALUE_VIZ_DATA "/seed_' + str(" FUNCTION_GET_SEED "()).zfill(5) + '/Scene";
const char* const DEFAULT_RXN_OUTPUT_FILENAME_PREFIX = "./" VALUE_REACT_DATA "/seed_' + str(" FUNCTION_GET_SEED "()).zfill(5) + '/";

const char* const INTERPRETER = "#!/usr/bin/env python3\n\n";

const char* const BASE_MODEL_IMPORTS =
    "import sys\n"
    "import os\n\n"
;

const char* const MCELL_DIR_SETUP =
    "MCELL_DIR = os.environ.get('MCELL_DIR', '')\n"
    "if MCELL_DIR:\n"
    "    sys.path.append(os.path.join(MCELL_DIR, 'lib'))\n"
    "else:\n"
    "    print(\"Error: variable MCELL_DIR that is used to find the mcell library was not set.\")\n"
    "    sys.exit(1)\n\n"
;

const char* const MCELL_IMPORT = "import mcell as m\n\n";

const char* const REGION_ALL_NAME = "ALL";
const char* const REGION_ALL_SUFFIX = "[ALL]";
const char* const NULL_PRODUCTS = "NULL";

const char* const COUNT = "COUNT";
const char* const WORLD = "WORLD";
const char* const WORLD_FIRST_UPPER = "World";

const char* const REV_RXN_SUFFIX = "_rev";
} // namespace MCell

#endif // SRC4_PYMCELLCONVERTER_GENERATOR_CONSTANTSS_H_
