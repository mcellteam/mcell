"""
Copyright (C) 2020 by
The Salk Institute for Biological Studies 

Use of this source code is governed by an MIT-style
license that can be found in the LICENSE file or at
https://opensource.org/licenses/MIT.
"""

# categories, defines input file names and 
# also ordering in documentation
CATEGORY_CONSTANTS = 'constants' 

CATEGORIES = [
    CATEGORY_CONSTANTS,
    'model',
    'simulation_setup',
    'subsystem', 
    'geometry',
    'instantiation',
    'observables',
    'callbacks',
    'introspection',
    'checkpointing',
    'submodules',
]


COPYRIGHT = \
"""/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/\n
"""

EXAMPLES_BASE_URL = 'https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/'

TARGET_DIRECTORY = '..' + '/' + 'generated'
DOC_DIRECTORY = '..'  + '/' + 'documentation'  + '/' + 'generated'
API_DIRECTORY = 'api'
WORK_DIRECTORY = '..'  + '/' + 'work'

KEY_ITEMS = 'items'
KEY_NAME = 'name'
ATTR_NAME_NAME = 'name' # attrribute with name 'name' is already defined in BaseDataClass
KEY_TYPE = 'type'
KEY_VALUE = 'value'
KEY_VALUES = 'values'
KEY_DEFAULT = 'default'
KEY_DOC = 'doc'
KEY_EXAMPLES = 'examples'
KEY_CATEGORY = 'category'  # set from input file name

KEY_SUPERCLASS = 'superclass'
KEY_SUPERCLASSES = 'superclasses'

KEY_CONSTANTS = 'constants'
KEY_ENUMS = 'enums'

KEY_METHODS = 'methods'
KEY_PARAMS = 'params'
KEY_RETURN_TYPE = 'return_type'
KEY_IS_CONST = 'is_const'

KEY_INHERITED = 'inherited' # used only internally, not in input YAML
KEY_NOT_AS_CTOR_ARG = 'not_as_ctor_arg'

VALUE_CLASS = 'class'
VALUE_SUBMODULE = 'submodule'

YAML_TYPE_PY_OBJECT = 'py::object'
YAML_TYPE_FLOAT = 'float'
YAML_TYPE_STR = 'str'
YAML_TYPE_INT = 'int'
YAML_TYPE_UINT32 = 'uint32'
YAML_TYPE_UINT64 = 'uint64'
YAML_TYPE_BOOL = 'bool'
YAML_TYPE_VEC2 = 'Vec2'
YAML_TYPE_VEC3 = 'Vec3'
YAML_TYPE_IVEC3 = 'IVec3'
YAML_TYPE_LIST = 'List'
YAML_TYPE_DICT = 'Dict'
YAML_TYPE_SPECIES = 'Species'
YAML_TYPE_ORIENTATION = 'Orientation'
YAML_TYPE_FUNCTION = 'std::function'

PYBIND_TYPE_OBJECT = 'object'
PY_CAST = 'py::cast'

CPP_TYPE_FLOAT = 'float'
CPP_TYPE_DOUBLE = 'double'
CPP_TYPE_STR = 'std::string'
CPP_TYPE_INT = 'int'
CPP_TYPE_UINT32 = 'uint'
CPP_TYPE_UINT64 = 'uint64_t'
CPP_TYPE_BOOL = 'bool'
CPP_TYPE_VEC2 = 'Vec2'
CPP_TYPE_VEC3 = 'Vec3'
CPP_TYPE_IVEC3 = 'IVec3'
CPP_VECTOR_TYPE = 'std::vector'
CPP_MAP_TYPE = 'std::map'

CPP_NONREFERENCE_TYPES = [CPP_TYPE_DOUBLE, CPP_TYPE_INT, CPP_TYPE_UINT64, CPP_TYPE_UINT32, CPP_TYPE_BOOL] # + CPP_VECTOR_TYPE but it must be checked with .startswith
CPP_REFERENCE_TYPES = [CPP_TYPE_STR, CPP_TYPE_VEC2, CPP_TYPE_VEC3, CPP_TYPE_IVEC3, CPP_VECTOR_TYPE]

UNSET_VALUE = 'unset'
EMPTY_ARRAY = 'empty'

UNSET_VALUE_FLOAT = 'FLT_UNSET'
UNSET_VALUE_STR = 'STR_UNSET'
UNSET_VALUE_INT = 'INT_UNSET'
UNSET_VALUE_UINT64 = '0' # default value, not unset
UNSET_VALUE_UINT32 = '0' # default value, not unset
UNSET_VALUE_VEC2 = 'VEC2_UNSET'
UNSET_VALUE_VEC3 = 'VEC3_UNSET'
UNSET_VALUE_ORIENTATION = 'Orientation::NOT_SET'
UNSET_VALUE_PTR = 'nullptr'

PY_NONE = 'None'

GEN_PREFIX = 'gen_'
GEN_GUARD_PREFIX = 'API_GEN_'
API_GUARD_PREFIX = 'API_'
GUARD_SUFFIX = '_H'
CTOR_SUFFIX = '_CTOR'
CTOR_NOARGS_SUFFIX = '_CTOR_NOARGS'
EXT_CPP = 'cpp'
EXT_H = 'h'

GEN_CLASS_PREFIX = 'Gen'
BASE_DATA_CLASS = 'BaseDataClass'
BASE_INTROSPECTION_CLASS = 'BaseIntrospectionClass'
BASE_EXPORT_CLASS = 'BaseExportClass'
CLASS_NAME_MODEL = 'Model'
CLASS_NAME_SUBSYSTEM = 'Subsystem'
CLASS_NAME_INSTANTIATION = 'Instantiation'
CLASS_NAME_OBSERVABLES = 'Observables'

CTOR_POSTPROCESS = 'postprocess_in_ctor'
RET_CTOR_POSTPROCESS = 'void' 

COPY_NAME = 'copy'
DEEPCOPY_NAME = 'deepcopy'
DEEPCOPY_VEC = 'deepcopy_vec'
DEEPCOPY_VEC_VEC = 'deepcopy_vec_vec'

IS_SET = 'is_set'

DEFAULT_CTOR_ARG_TYPE = 'DefaultCtorArgType'

RET_TYPE_CHECK_SEMANTICS = 'void'
CHECK_SEMANTICS = 'check_semantics'
DECL_CHECK_SEMANTICS = CHECK_SEMANTICS + '() const'
DECL_DEFINE_PYBINDIND_CONSTANTS = 'void define_pybinding_constants(py::module& m)'
DECL_DEFINE_PYBINDIND_ENUMS = 'void define_pybinding_enums(py::module& m)'
DECL_SET_INITIALIZED = 'set_initialized()'

RET_TYPE_SET_ALL_DEFAULT_OR_UNSET = 'void'
SET_ALL_DEFAULT_OR_UNSET_DECL = 'set_all_attributes_as_default_or_unset()'
SET_ALL_CUSTOM_TO_DEFAULT_DECL = 'set_all_custom_attributes_to_default()'

RET_TYPE_TO_STR = 'std::string'
SHARED_PTR = 'std::shared_ptr'
MAKE_SHARED = 'std::make_shared'
DECL_TO_STR_W_DEFAULT = 'to_str(const bool all_details=false, const std::string ind="") const'
DECL_TO_STR = 'to_str(const bool all_details, const std::string ind) const'

RET_TYPE_EXPORT_TO_PYTHON = 'std::string' 

CTX = 'ctx'
EXPORTED_NAME = 'exported_name'
EXPORT_TO_PYTHON_ARGS = 'std::ostream& out, PythonExportContext& ' + CTX
DECL_EXPORT_TO_PYTHON = 'export_to_python(' + EXPORT_TO_PYTHON_ARGS + ')'
EXPORT_VEC_PREFIX = 'export_vec_'
M_DOT = 'm.'

KEYWORD_OVERRIDE = 'override'
KEYWORD_VIRTUAL = 'virtual'
  
CLASS_NAME_ATTR = 'class_name'
CACHED_DATA_ARE_UPTODATE_ATTR = 'cached_data_are_uptodate'

GEN_VECTORS_BIND = 'gen_vectors_bind(py::module& m)'

GEN_CONSTANTS_H = 'gen_constants.h'
GEN_CONSTANTS_CPP = 'gen_constants.cpp'

GEN_VECTORS_MAKE_OPAQUE_H = 'gen_vectors_make_opaque.h'
GEN_VECTORS_BIND_CPP = 'gen_vectors_bind.cpp'

MCELL_PYI = 'mcell.pyi'
EXT_RST = '.rst'
API_RST = 'api' + EXT_RST 

GEN_NAMES_H = 'gen_names.h'
NAME_PREFIX = 'NAME_'
CLASS_PREFIX = 'CLASS_'
ENUM_PREFIX = 'ENUM_'
ENUM_VALUE_PREFIX = 'EV_'
CONSTANT_VALUE_PREFIX = 'CV_'

INCLUDE_API_MCELL_H = '#include "api/mcell.h"'
INCLUDE_API_COMMON_H = '#include "api/api_common.h"'
INCLUDE_API_PYTHON_EXPORT_UTILS_H = '#include "api/python_export_utils.h"'
INCLUDE_API_BASE_DATA_CLASS_H = '#include "api/base_data_class.h"'
INCLUDE_API_BASE_EXPORT_CLASS_H = '#include "api/base_export_class.h"'
INCLUDE_API_BASE_INTROSPECTION_CLASS_H = '#include "api/base_introspection_class.h"'
NAMESPACES_BEGIN = 'namespace MCell {\nnamespace API {'
NAMESPACES_END = '} // namespace API\n} // namespace MCell'

VEC_NONPTR_TO_STR = 'vec_nonptr_to_str'
VEC_PTR_TO_STR = 'vec_ptr_to_str'
F_TO_STR = 'f_to_str'

PY_BIND_VECTOR = 'py::bind_vector'
PY_IMPLICITLY_CONVERTIBLE = 'py::implicitly_convertible'

PYBIND11_MAKE_OPAQUE = 'PYBIND11_MAKE_OPAQUE'
