#ifndef _DATAMODEL_DEFINES_H_
#define _DATAMODEL_DEFINES_H_

#include "defines.h"
#include "json/json.h"

const char* const DEFAULT_DATAMODEL_FILENAME = "datamodel.json";

// used also outside of datamodel, maybe move elsewhere
const char* const ALL_MOLECULES = "ALL_MOLECULES";
const char* const ALL_VOLUME_MOLECULES = "ALL_VOLUME_MOLECULES";
const char* const ALL_SURFACE_MOLECULES = "ALL_SURFACE_MOLECULES";

// ---------------------------------- datamodel constants----------------------------------

const char* const JSON_DM_VERSION_1638 = "DM_2014_10_24_1638";
const char* const JSON_DM_VERSION_1632 = "DM_2018_10_16_1632";
const char* const JSON_DM_VERSION_1330 = "DM_2018_01_11_1330";
const char* const JSON_DM_VERSION_1700 = "DM_2015_04_13_1700";
const char* const JSON_DM_VERSION_1300 = "DM_2017_06_23_1300";

const int BLENDER_VERSION[] = {2, 79, 0};

const char* const KEY_ROOT = "root"; // not used, mainly for error reporting

// all keys should be defined as a constant string
const char* const KEY_MCELL = "mcell";

const char* const KEY_NAME = "name";
const char* const VALUE_NONE = "None";
const char* const VALUE_BLENDER = "blender";

const char* const KEY_DATA_MODEL_VERSION = "data_model_version";
const char* const KEY_CELLBLENDER_VERSION = "cellblender_version";
const char* const VALUE_CELLBLENDER_VERSION = "0.1.54";


const char* const KEY_GEOMETRICAL_OBJECTS = "geometrical_objects";
const char* const KEY_OBJECT_LIST = "object_list";
const char* const KEY_LOCATION = "location";
const char* const KEY_ELEMENT_CONNECTIONS = "element_connections";
const char* const KEY_VERTEX_LIST = "vertex_list";
const char* const KEY_DEFINE_SURFACE_REGIONS = "define_surface_regions";
const char* const KEY_INCLUDE_ELEMENTS = "include_elements";

const char* const KEY_REACTION_LIST = "reaction_list";
const char* const KEY_DEFINE_REACTIONS = "define_reactions";
const char* const KEY_RXN_NAME = "rxn_name";
const char* const KEY_RXN_TYPE = "rxn_type";
const char* const VALUE_REVERSIBLE = "reversible";
const char* const VALUE_IRREVERSIBLE = "irreversible";
const char* const KEY_REACTANTS = "reactants";
const char* const KEY_PRODUCTS = "products";
const char* const KEY_FWD_RATE= "fwd_rate";
const char* const KEY_BKWD_RATE= "bkwd_rate";
const char* const KEY_VARIABLE_RATE_VALID = "variable_rate_valid";
const char* const KEY_VARIABLE_RATE = "variable_rate";
const char* const KEY_VARIABLE_RATE_SWITCH = "variable_rate_switch";
const char* const KEY_VARIABLE_RATE_TEXT = "variable_rate_text";

const char* const KEY_MATERIAL_NAMES = "material_names";
const char* const KEY_MATERIALS = "materials";
const char* const KEY_MATERIAL_DICT = "material_dict";
const char* const KEY_VALUE_MEMBRANE = "membrane"; // used both as a key and value
const char* const KEY_DIFFUSE_COLOR = "diffuse_color";
const char* const KEY_R = "r";
const char* const KEY_G = "g";
const char* const KEY_B = "b";
const char* const KEY_A = "a";

const char* const KEY_MODEL_OBJECTS = "model_objects";
const char* const KEY_MODEL_OBJECT_LIST = "model_object_list";
const char* const KEY_PARENT_OBJECT = "parent_object";
const char* const KEY_DESCRIPTION = "description";
const char* const KEY_OBJECT_SOURCE = "object_source";
const char* const KEY_DYNAMIC_DISPLAY_SOURCE = "dynamic_display_source";
const char* const KEY_SCRIPT_NAME = "script_name";
const char* const KEY_MEMBRANE_NAME = "membrane_name";
const char* const KEY_DYNAMIC = "dynamic";

const char* const KEY_DEFINE_SURFACE_CLASSES = "define_surface_classes";

const char* const KEY_MOL_VIZ = "mol_viz";
const char* const KEY_MANUAL_SELECT_VIZ_DIR = "manual_select_viz_dir";
const char* const KEY_FILE_START_INDEX = "file_start_index";
const char* const KEY_SEED_LIST = "seed_list";
const char* const KEY_COLOR_LIST = "color_list";
const char* const KEY_ACTIVE_SEED_INDEX = "active_seed_index";
const char* const KEY_FILE_INDEX = "file_index";
const char* const KEY_FILE_NUM = "file_num";
const char* const KEY_VIZ_ENABLE = "viz_enable";
const char* const KEY_FILE_NAME = "file_name";
const char* const KEY_COLOR_INDEX = "color_index";
const char* const KEY_RENDER_AND_SAVE = "render_and_save";
const char* const KEY_FILE_STEP_INDEX = "file_step_index";
const char* const KEY_VIZ_LIST = "viz_list";
const char* const KEY_FILE_STOP_INDEX = "file_stop_index";
const char* const KEY_FILE_DIR = "file_dir";

const char* const KEY_MODEL_LANGUAGE = "model_language";
const char* const VALUE_MCELL3 = "mcell3";
const char* const KEY_PARAMETER_SYSTEM = "parameter_system";

const char* const KEY_INITIALIZATION = "initialization";

const char* const KEY_BLENDER_VERSION = "blender_version";

const char* const KEY_MOLECULE_LIST = "molecule_list";
const char* const KEY_DEFINE_MOLECULES = "define_molecules";
const char* const KEY_DISPLAY = "display";
const char* const KEY_EMIT = "emit";
const char* const KEY_COLOR = "color";
const char* const KEY_GLYPH = "glyph";
const char* const VALUE_GLYPH_SPHERE_1 = "Sphere_1";
const char* const KEY_SCALE = "scale";
const char* const KEY_BNGL_COMPONENT_LIST = "bngl_component_list";
const char* const KEY_MOL_BNGL_LABEL = "mol_bngl_label";
const char* const KEY_MOL_NAME = "mol_name";
const char* const KEY_MOL_TYPE = "mol_type";
const char* const VALUE_MOL_TYPE_2D = "2D";
const char* const VALUE_MOL_TYPE_3D = "3D";
const char* const KEY_DIFFUSION_CONSTANT = "diffusion_constant";
const char* const KEY_EXPORT_VIZ = "export_viz";
const char* const KEY_CUSTOM_SPACE_STEP = "custom_space_step";
const char* const KEY_MAXIMUM_STEP_LENGTH = "maximum_step_length";
const char* const KEY_TARGET_ONLY = "target_only";
const char* const KEY_SPATIAL_STRUCTURE = "spatial_structure";
const char* const KEY_CUSTOM_TIME_STEP = "custom_time_step";

const char* const KEY_RELEASE_SITES = "release_sites";
const char* const KEY_RELEASE_SITE_LIST = "release_site_list";
const char* const KEY_POINTS_LIST = "points_list";
const char* const KEY_PATTERN = "pattern";
const char* const KEY_STDDEV = "stddev";
const char* const KEY_QUANTITY = "quantity";
const char* const KEY_MOLECULE = "molecule";
const char* const KEY_SHAPE = "shape";
const char* const VALUE_OBJECT = "OBJECT";
const char* const KEY_RELEASE_PROBABILITY = "release_probability";
const char* const KEY_OBJECT_EXPR = "object_expr";
const char* const KEY_ORIENT = "orient";
const char* const KEY_LOCATION_X = "location_x";
const char* const KEY_LOCATION_Y = "location_y";
const char* const KEY_LOCATION_Z = "location_z";
const char* const KEY_SITE_DIAMETER = "site_diameter";
const char* const KEY_QUANTITY_TYPE = "quantity_type";
const char* const VALUE_NUMBER_TO_RELEASE = "NUMBER_TO_RELEASE";
const char* const VALUE_GAUSSIAN_RELEASE_NUMBER = "GAUSSIAN_RELEASE_NUMBER";
const char* const VALUE_DENSITY = "DENSITY";


const double DEFAULT_OBJECT_COLOR_COMPONENT = 0.8;
const double DEFAULT_OBJECT_ALPHA = 0.25;


// ---------------------------------- datamodel utilities----------------------------------

namespace DMUtil {

static inline void json_add_version(Json::Value& define_molecules, const char* ver) {
  define_molecules[KEY_DATA_MODEL_VERSION] = ver;
}

static inline void json_append_triplet(Json::Value& list, const float x, const float y, const float z) {
  Json::Value list_triplet;
  list_triplet.append(Json::Value(x));
  list_triplet.append(Json::Value(y));
  list_triplet.append(Json::Value(z));
  list.append(list_triplet);
}

static inline std::string remove_obj_name_prefix(const std::string& prefix, const std::string& name) {
  size_t pos = prefix.size() + 1;
  assert(name.size() > pos);
  assert(name.substr(0, pos) == prefix + ".");
  return name.substr(pos);
}

// we do not know the prefix in this case
static inline std::string remove_obj_name_prefix(const std::string& name) {
  size_t pos = name.find('.');
  assert(pos != std::string::npos);
  assert(pos + 1 < name.size());
  return name.substr(pos + 1);
}

static inline std::string get_object_w_region_name(const std::string& name) {

  std::string noprefix = remove_obj_name_prefix(name);
  size_t pos = noprefix.find(',');
  assert(pos != std::string::npos);
  assert(pos + 1 < name.size());

  std::string obj = noprefix.substr(0, pos);
  std::string reg = noprefix.substr(pos + 1);
  return obj + "[" + reg + "]";
}

static inline std::string get_surface_region_name(const std::string& name) {
  size_t pos = name.find(',');
  // TODO: these checks should be enabled also for release build
  assert(pos != std::string::npos);
  assert(pos + 1 < name.size());
  return name.substr(pos + 1);
}

static inline std::string orientation_to_str(const orientation_t o) {
  switch (o) {
    case ORIENTATION_DOWN: return ",";
    case ORIENTATION_UP: return "'";
    case ORIENTATION_NONE: return ";"; // this character is used by cellblender for volume mols is release sites
    default:
      assert(false);
      return "ERROR";
  }
}

#define CONVERSION_UNSUPPORTED(msg) \
  do { errs() << msg << " This is not supported yet.\n"; exit(1); } while (0)

#define CONVERSION_CHECK(cond, msg) \
  do { if (!(cond)) { errs() << "Check failed: " << #cond << ": " << msg << " This is not supported yet.\n"; exit(1); } } while (0)

} // namespace DataModel

#endif // _DATAMODEL_DEFINES_H_
