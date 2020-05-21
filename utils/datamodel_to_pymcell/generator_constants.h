
#ifndef SRC4_PYMCELLCONVERTER_GENERATOR_CONSTANTSS_H_
#define SRC4_PYMCELLCONVERTER_GENERATOR_CONSTANTSS_H_

namespace MCell {

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

const char* const PARAM_SEED = "SEED";
const char* const PARAM_ITERATIONS = "ITERATIONS";
const char* const PARAM_TIME_STEP = "TIME_STEP";
const char* const PARAM_DUMP = "DUMP";

const char* const VEC3 = "Vec3";

const char* const UNNAMED_REACTION_RULE_PREFIX = "unnamed_reaction_rule_";

const char* const VIZ_OUTPUT_NAME = "viz_output";
const char* const DEFAULT_VIZ_OUTPUT_FILENAME_PREFIX = "'./viz_data/seed_' + str(SEED).zfill(5) + '/Scene'";

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

} // namespace MCell

#endif // SRC4_PYMCELLCONVERTER_GENERATOR_CONSTANTSS_H_
