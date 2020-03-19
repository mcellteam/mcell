#ifndef _DATAMODEL_DEFINES_H_
#define _DATAMODEL_DEFINES_H_

#include "json/json.h"

// ---------------------------------- datamodel constants----------------------------------

const char* const JSON_DM_VERSION_1638 = "DM_2014_10_24_1638";
const char* const JSON_DM_VERSION_1632 = "DM_2018_10_16_1632";
const char* const JSON_DM_VERSION_1330 = "DM_2018_01_11_1330";
const char* const JSON_DM_VERSION_1700 = "DM_2015_04_13_1700";
const char* const JSON_DM_VERSION_1300 = "DM_2017_06_23_1300";

const int BLENDER_VERSION[] = {2, 79, 0};


// TODO: move this into a different location
const char* const ALL_MOLECULES = "ALL_MOLECULES";
const char* const ALL_VOLUME_MOLECULES = "ALL_VOLUME_MOLECULES";
const char* const ALL_SURFACE_MOLECULES = "ALL_SURFACE_MOLECULES";


// all keys should be defined as a constant string
const char* const KEY_MCELL = "mcell";

const char* const KEY_NAME = "name";
const char* const VALUE_NONE = "None";
const char* const VALUE_BLENDER = "blender";

const char* const KEY_DATA_MODEL_VERSION = "data_model_version";
const char* const KEY_CELLBLENDER_VERSION = "cellblender_version";
const char* const VALUE_CELLBLENDER_VERSION = "0.1.54";

const char* const KEY_MOLECULE_LIST = "molecule_list";

const char* const KEY_GEOMETRICAL_OBJECTS = "geometrical_objects";
const char* const KEY_OBJECT_LIST = "object_list";
const char* const KEY_ELEMENT_CONNECTIONS = "element_connections";
const char* const KEY_VERTEX_LIST = "vertex_list";

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

const double DEFAULT_OBJECT_COLOR_COMPONENT = 0.8;
const double DEFAULT_OBJECT_ALPHA = 0.25;


// ---------------------------------- datamodel utilities----------------------------------

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

#endif // _DATAMODEL_DEFINES_H_
