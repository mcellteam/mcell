#!/usr/bin/env python3
"""
Copyright (C) 2020 by
The Salk Institute for Biological Studies 

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

For the complete terms of the GNU General Public License, please see this URL:
http://www.gnu.org/licenses/gpl-2.0.html
"""

# TODO: change 'items' to 'attributes'?
# TODO: split to at least two files
# TODO: cleanup const/nonconst handling
# TODO: unify superclass vs superclasses - 
#       w superclasses, the objects are inherited manually now and it should be made automatic
#       or maybe just rename it..
# TODO: print unset constants as UNSET, not diffusion_constant_2d=3.40282e+38

import sys
import os
import yaml
import re
from datetime import datetime
from copy import copy

VERBOSE = False # may be overridden by argument -v 

# order is not important
ALL_INPUT_FILES = [
    'constants.yaml',
    'subsystem.yaml', 
    'geometry.yaml',
    'instantiation.yaml',
    'observables.yaml',
    'simulation_setup.yaml',
    'callbacks.yaml',
    'model.yaml',
    'introspection.yaml',
    'submodules.yaml'
]


COPYRIGHT = \
"""/******************************************************************************
 *
 * Copyright (C) 2020 by
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
******************************************************************************/\n
"""

TARGET_DIRECTORY = os.path.join('..', 'generated')
DOC_DIRECTORY = os.path.join('..', 'documentation', 'generated')
API_DIRECTORY = os.path.join('..', 'api')
WORK_DIRECTORY = os.path.join('..', 'work')

KEY_ITEMS = 'items'
KEY_NAME = 'name'
ATTR_NAME_NAME = 'name' # attrribute with name 'name' is already defined in BaseDataClass
KEY_TYPE = 'type'
KEY_VALUE = 'value'
KEY_VALUES = 'values'
KEY_DEFAULT = 'default'
KEY_DOC = 'doc'

KEY_SUPERCLASS = 'superclass'
KEY_SUPERCLASSES = 'superclasses'

KEY_CONSTANTS = 'constants'
KEY_ENUMS = 'enums'

KEY_METHODS = 'methods'
KEY_PARAMS = 'params'
KEY_RETURN_TYPE = 'return_type'

KEY_INHERITED = 'inherited' # used only internally, not in input YAML
KEY_NOT_AS_CTOR_ARG = 'not_as_ctor_arg'

VALUE_CLASS = 'class'
VALUE_SUBMODULE = 'submodule'

YAML_TYPE_PY_OBJECT = 'py::object'
YAML_TYPE_FLOAT = 'float'
YAML_TYPE_STR = 'str'
YAML_TYPE_INT = 'int'
YAML_TYPE_LONG = 'long'
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

CPP_TYPE_FLOAT = 'float_t'
CPP_TYPE_STR = 'std::string'
CPP_TYPE_INT = 'int'
CPP_TYPE_LONG = 'long'
CPP_TYPE_BOOL = 'bool'
CPP_TYPE_VEC2 = 'Vec2'
CPP_TYPE_VEC3 = 'Vec3'
CPP_TYPE_IVEC3 = 'IVec3'
CPP_VECTOR_TYPE = 'std::vector'
CPP_MAP_TYPE = 'std::map'

CPP_NONREFERENCE_TYPES = [CPP_TYPE_FLOAT, CPP_TYPE_INT, CPP_TYPE_LONG, CPP_TYPE_BOOL] # + CPP_VECTOR_TYPE but it must be checked with .startswith
CPP_REFERENCE_TYPES = [CPP_TYPE_STR, CPP_TYPE_VEC2, CPP_TYPE_VEC3, CPP_TYPE_IVEC3, CPP_VECTOR_TYPE]

UNSET_VALUE = 'unset'
EMPTY_ARRAY = 'empty'

UNSET_VALUE_FLOAT = 'FLT_UNSET'
UNSET_VALUE_STR = 'STR_UNSET'
UNSET_VALUE_INT = 'INT_UNSET'
UNSET_VALUE_LONG = 'LONG_UNSET'
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

CTOR_POSTPROCESS = 'postprocess_in_ctor'
RET_CTOR_POSTPROCESS = 'void' 

RET_TYPE_CHECK_SEMANTICS = 'void'
CHECK_SEMANTICS = 'check_semantics'
DECL_CHECK_SEMANTICS = CHECK_SEMANTICS + '() const'
DECL_DEFINE_PYBINDIND_CONSTANTS = 'void define_pybinding_constants(py::module& m)'
DECL_DEFINE_PYBINDIND_ENUMS = 'void define_pybinding_enums(py::module& m)'
DECL_SET_INITIALIZED = 'set_initialized()'

RET_TYPE_SET_ALL_DEFAULT_OR_UNSET = 'void'
SET_ALL_DEFAULT_OR_UNSET_DECL = 'set_all_attributes_as_default_or_unset()'

RET_TYPE_TO_STR = 'std::string'
SHARED_PTR = 'std::shared_ptr'
DECL_TO_STR_W_DEFAULT = 'to_str(const std::string ind="") const'
DECL_TO_STR = 'to_str(const std::string ind) const'
KEYWORD_OVERRIDE = 'override'
  
CLASS_NAME_ATTR = 'class_name'

GEN_CONSTANTS_H = 'gen_constants.h'
GEN_CONSTANTS_CPP = 'gen_constants.cpp'

MCELL_PYI = 'mcell.pyi'
CLASSES_RST = 'classes.rst'

GEN_NAMES_H = 'gen_names.h'
NAME_PREFIX = 'NAME_'
CLASS_PREFIX = 'CLASS_'
ENUM_PREFIX = 'ENUM_'
ENUM_VALUE_PREFIX = 'EV_'
CONSTANT_VALUE_PREFIX = 'CV_'

INCLUDE_API_MCELL_H = '#include "../api/mcell.h"'
INCLUDE_API_COMMON_H = '#include "../api/common.h"'
INCLUDE_API_BASE_DATA_CLASS_H = '#include "../api/base_data_class.h"'
INCLUDE_API_BASE_INTROSPECTION_CLASS_H = '#include "../api/base_introspection_class.h"'
NAMESPACES_BEGIN = 'namespace MCell {\nnamespace API {'
NAMESPACES_END = '} // namespace API\n} // namespace MCell'

VEC_NONPTR_TO_STR = 'vec_nonptr_to_str'
VEC_PTR_TO_STR = 'vec_ptr_to_str'

# --- some global data ---
g_enums = set()

# ------------------------

def rename_cpp_reserved_id(name):
    if (name == 'union'):
        return 'union_'
    else:
        return name


def get_underscored(class_name):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', class_name).lower()

def get_gen_header_guard_name(class_name):
    return GEN_GUARD_PREFIX + get_underscored(class_name).upper() + GUARD_SUFFIX

def get_api_header_guard_name(class_name):
    return API_GUARD_PREFIX + get_underscored(class_name).upper() + GUARD_SUFFIX

def get_gen_class_file_name(class_name, extension):
    return os.path.join(GEN_PREFIX + get_underscored(class_name) + '.' + extension)

def get_gen_class_file_name_w_dir(class_name, extension):
    return os.path.join(TARGET_DIRECTORY, get_gen_class_file_name(class_name, extension))

def get_api_class_file_name(class_name, extension):
    return get_underscored(class_name) + '.' + extension

def get_api_class_file_name_w_dir(class_name, extension):
    return os.path.join(API_DIRECTORY, get_api_class_file_name(class_name, extension))
   
def get_api_class_file_name_w_work_dir(class_name, extension):
    return os.path.join(WORK_DIRECTORY, get_api_class_file_name(class_name, extension))

def is_yaml_list_type(t):
    return t.startswith(YAML_TYPE_LIST)

def is_yaml_dict_type(t):
    return t.startswith(YAML_TYPE_DICT)

def is_yaml_function_type(t):
    return t.startswith(YAML_TYPE_FUNCTION)

# rename inner to underlying?
def get_inner_list_type(t):
    if is_yaml_list_type(t):
        return t[len(YAML_TYPE_LIST)+1:-1]
    else:
        return t

def get_inner_dict_key_type(t):
    if is_yaml_dict_type(t):
        return t.replace('[', ',').replace(']', ',').split(',')[1].strip()
    else:
        return t

def get_inner_dict_value_type(t):
    if is_yaml_dict_type(t):
        return t.replace('[', ',').replace(']', ',').split(',')[2].strip()
    else:
        return t
    
def get_inner_function_type(t):
    # callbacks pass back only one shared_ptr argument 
    if is_yaml_function_type(t):
        return t.split('<')[2].split('>')[0].strip()
    else:
        return t    

def get_first_inner_type(t):
    if is_yaml_list_type(t):
        return get_inner_list_type(t)
    elif is_yaml_dict_type(t):
        return get_inner_dict_key_type(t)
    elif is_yaml_function_type(t):
        return get_inner_function_type(t)
    else: 
        return t    

def is_base_yaml_type(t):
    return \
        t == YAML_TYPE_FLOAT or t == YAML_TYPE_STR or t == YAML_TYPE_INT or t == YAML_TYPE_LONG or \
        t == YAML_TYPE_BOOL or t == YAML_TYPE_VEC2 or t == YAML_TYPE_VEC3 or t == YAML_TYPE_IVEC3 or \
        t == YAML_TYPE_PY_OBJECT or \
        (is_yaml_function_type(t) and is_base_yaml_type(get_inner_function_type(t))) or \
        (is_yaml_list_type(t) and is_base_yaml_type(get_inner_list_type(t))) or \
        (is_yaml_dict_type(t) and is_base_yaml_type(get_inner_dict_key_type(t)) and is_base_yaml_type(get_inner_dict_value_type(t)))


def is_yaml_ptr_type(t):
    if t == YAML_TYPE_PY_OBJECT:
        return False
    else:
        return t[-1] == '*'


def yaml_type_to_cpp_type(t):
    assert len(t) >= 1
    if t == YAML_TYPE_FLOAT:
        return CPP_TYPE_FLOAT
    elif t == YAML_TYPE_STR:
        return CPP_TYPE_STR
    elif t == YAML_TYPE_INT:
        return CPP_TYPE_INT
    elif t == YAML_TYPE_LONG:
        return CPP_TYPE_LONG
    elif t == YAML_TYPE_LONG:
        return CPP_TYPE_BOOL
    elif t == YAML_TYPE_VEC2:
        return CPP_TYPE_VEC2
    elif t == YAML_TYPE_VEC3:
        return CPP_TYPE_VEC3
    elif t == YAML_TYPE_IVEC3:
        return CPP_TYPE_IVEC3
    elif is_yaml_list_type(t):
        assert len(t) > 7
        inner_type = yaml_type_to_cpp_type(get_inner_list_type(t))
        return CPP_VECTOR_TYPE + '<' + inner_type + '>'
    elif is_yaml_dict_type(t):
        key_type = yaml_type_to_cpp_type(get_inner_dict_key_type(t))
        value_type = yaml_type_to_cpp_type(get_inner_dict_value_type(t))
        return CPP_MAP_TYPE + '<' + key_type + ', ' + value_type + '>'
        
    else:
        if is_yaml_ptr_type(t):
            return SHARED_PTR + '<' + t[0:-1] + '>' 
        else:
            return t # standard ttype
    
    
def get_cpp_bool_string(val):
    if val == 'True':
        return 'true'
    elif val == 'False':
        return 'false'
    else:
        assert false

    
def yaml_type_to_pybind_type(t):
    assert len(t) >= 1
    if t == YAML_TYPE_FLOAT:
        return 'float_' 
    elif t == YAML_TYPE_STR:
        return YAML_TYPE_STR
    elif t == YAML_TYPE_BOOL:
        return CPP_TYPE_BOOL
    elif t == YAML_TYPE_INT:
        return CPP_TYPE_INT + '_'
    elif t == YAML_TYPE_LONG:
        return CPP_TYPE_INT + '_'
    elif t == YAML_TYPE_SPECIES:
        return PYBIND_TYPE_OBJECT
    else:
        assert False, "Unsupported constant type " + t 
    
    
def yaml_type_to_py_type(t):
    assert len(t) >= 1
    if t == YAML_TYPE_LONG:
        return YAML_TYPE_INT # not sure what should be the name
    elif is_yaml_function_type(t):
        return "Callable, # " + t
    elif t == YAML_TYPE_PY_OBJECT:
        return "Any, # " + t
    else:
        return t.replace('*', '')

        
def is_cpp_ptr_type(cpp_type):
    return cpp_type.startswith(SHARED_PTR)

    
def is_cpp_ref_type(cpp_type):
    not_reference = \
            cpp_type in CPP_NONREFERENCE_TYPES or \
            is_cpp_ptr_type(cpp_type) or \
            cpp_type.startswith(CPP_VECTOR_TYPE) or \
            cpp_type in g_enums or \
            is_yaml_function_type(cpp_type) or \
            cpp_type == YAML_TYPE_PY_OBJECT
    
    return not not_reference
    
    
def get_type_as_ref_param(attr):
    assert KEY_TYPE in attr
    
    yaml_type = attr[KEY_TYPE]
    cpp_type = yaml_type_to_cpp_type(yaml_type)
    
    res = cpp_type
    if is_cpp_ref_type(cpp_type):
        res += '&'
    
    return res


def get_default_or_unset_value(attr):
    assert KEY_TYPE in attr
    t = attr[KEY_TYPE]

    if KEY_DEFAULT in attr:
        default_value = attr[KEY_DEFAULT]
        if default_value != UNSET_VALUE and default_value != EMPTY_ARRAY:
            res = str(default_value)
            # might need to convert enum.value into enum::value
            if not is_base_yaml_type(t):
                res = res.replace('.', '::')
            elif t == YAML_TYPE_BOOL:
                res = get_cpp_bool_string(res)
            elif t == YAML_TYPE_STR:
                # default strings must be quoted (not the 'unset' one because the UNSET_STR is a constant)
                res = '"' + res + '"'
                
            return res 
    
    if t == YAML_TYPE_FLOAT:
        return UNSET_VALUE_FLOAT
    elif t == YAML_TYPE_STR:
        return UNSET_VALUE_STR
    elif t == YAML_TYPE_INT:
        return UNSET_VALUE_INT
    elif t == YAML_TYPE_LONG:
        return UNSET_VALUE_LONG
    elif t == YAML_TYPE_BOOL:
        assert False, "There is no unset value for bool - for " + attr[KEY_NAME]
        return "error"
    elif t == YAML_TYPE_VEC2:
        return UNSET_VALUE_VEC2
    elif t == YAML_TYPE_VEC3:
        return UNSET_VALUE_VEC3
    elif t == YAML_TYPE_IVEC3:
        return UNSET_VALUE_IVEC3
    elif t == YAML_TYPE_ORIENTATION:
        return UNSET_VALUE_ORIENTATION
    elif is_yaml_list_type(t):
        return yaml_type_to_cpp_type(t) + '()'
    elif is_yaml_dict_type(t):
        return yaml_type_to_cpp_type(t) + '()'
    else:
        return UNSET_VALUE_PTR

def get_default_or_unset_value_py(attr):
    assert KEY_TYPE in attr
    t = attr[KEY_TYPE]

    if KEY_DEFAULT in attr:
        default_value = attr[KEY_DEFAULT]
        if default_value == "":
            return "''"
        elif default_value != UNSET_VALUE and default_value != EMPTY_ARRAY:
            return str(default_value)
            # might need to convert enum.value into enum::value
            #if not is_base_yaml_type(t):
            #    res = res.replace('.', '::')
            #elif t == YAML_TYPE_BOOL:
            #    res = get_cpp_bool_string(res)
    
    return PY_NONE


def write_generated_notice(f):
    now = datetime.now()
    # date printing is disabled during the development phase to minimize changes in files 
    #date_time = now.strftime("%m/%d/%Y, %H:%M")
    #f.write('// This file was generated automatically on ' + date_time + ' from ' + '\'' + input_file_name + '\'\n\n')

def write_ctor_decl(f, class_def, class_name, append_backslash, indent_and_fix_rst_chars, only_inherited, with_args=True):
    items = class_def[KEY_ITEMS] if KEY_ITEMS in class_def else []
    backshlash = '\\' if append_backslash else ''
    
    f.write(indent_and_fix_rst_chars + class_name + '( ' + backshlash + '\n')
     
    inherited_items = [ attr for attr in items if is_inherited(attr) ]

    if only_inherited:
        # override also the original items list
        items = inherited_items 
     
    # ctor parameters    
    if with_args:
        num_items = len(items)
        for i in range(num_items):
            attr = items[i]
            
            if not only_inherited and attr_not_in_ctor(attr):
                continue


            assert KEY_NAME in attr
            name = attr[KEY_NAME]
    
            const_spec = 'const ' if not is_yaml_ptr_type(attr[KEY_TYPE]) else ''
            f.write(indent_and_fix_rst_chars + '    ' + const_spec + get_type_as_ref_param(attr) + ' ' + name + '_')
            
            if KEY_DEFAULT in attr:
                f.write(' = ' + get_default_or_unset_value(attr))
            
            if i != num_items - 1:
                f.write(',')
            f.write(' ' + backshlash + '\n')
    
    f.write(indent_and_fix_rst_chars + ') ')
    
    if has_superclass_other_than_base(class_def):
        # call superclass ctor
        # only one, therefore all inherited attributes are its arguments
        if only_inherited:
            # we are generating ctor for the superclass of the Gen class, e.g. for GenSpecies 
            # and we need to initialize ComplexInstance 
            superclass_name = class_def[KEY_SUPERCLASS]
        else:
            # we are generating ctor for the superclass, e.g. for Species and we need to initialize 
            # GenSpecies
            superclass_name = GEN_CLASS_PREFIX + class_name
            
        f.write(' : ' + superclass_name + '(')
        
        num_inherited_items = len(inherited_items)
        for i in range(num_inherited_items):
            f.write(inherited_items[i][KEY_NAME] + '_')
            if i != num_inherited_items - 1:
                f.write(',')
                
        f.write(') ')
    

def write_ctor_define(f, class_def, class_name):
    with_args = True
    if KEY_SUPERCLASS in class_def and class_def[KEY_SUPERCLASS] == BASE_INTROSPECTION_CLASS:
        with_args = False
        
    suffix = CTOR_SUFFIX if with_args else CTOR_NOARGS_SUFFIX
    f.write('#define ' + get_underscored(class_name).upper() + suffix + '() \\\n')

    write_ctor_decl(f, class_def, class_name, append_backslash=True, indent_and_fix_rst_chars='    ', only_inherited=False, with_args=with_args)

    f.write('{ \\\n')

    # initialization code
    f.write('      ' + CLASS_NAME_ATTR + ' = "' + class_name + '"; \\\n')
    items = class_def[KEY_ITEMS] if KEY_ITEMS in class_def else []    
    num_items = len(items)
    for i in range(num_items):
        if attr_not_in_ctor(items[i]):
            continue

        assert KEY_NAME in items[i] 
        attr_name = items[i][KEY_NAME]
        if with_args:
            f.write('      ' + attr_name + ' = ' + attr_name + '_; \\\n')
        else:
            f.write('      ' + attr_name + ' = ' + get_default_or_unset_value(items[i]) + '; \\\n')
            
    f.write('      ' + CTOR_POSTPROCESS + '();\\\n')    
    f.write('      ' + CHECK_SEMANTICS + '();\\\n')
    f.write('    }\n\n')    
    
    
def write_ctor_for_superclass(f, class_def, class_name):
    write_ctor_decl(f, class_def, GEN_CLASS_PREFIX + class_name, append_backslash=False, indent_and_fix_rst_chars='  ', only_inherited=True)
    f.write(' {\n')
    f.write('  }\n')
    

def write_attr_with_get_set(f, class_def, attr):
    assert KEY_NAME in attr, KEY_NAME + " is not set in " + str(attr) 
    assert KEY_TYPE in attr, KEY_TYPE + " is not set in " + str(attr)
    
    name = attr[KEY_NAME]
    
    # skip attribute 'name'
    if name == ATTR_NAME_NAME:
        return False
    
    yaml_type = attr[KEY_TYPE]
    cpp_type = yaml_type_to_cpp_type(yaml_type)
    
    #decl_const = 'const ' if is_cpp_ptr_type(cpp_type) else ''
    decl_type = cpp_type
    
    # decl
    f.write('  ' + decl_type + ' ' + name + ';\n')
    
    # setter
    arg_type_const = 'const ' if not is_cpp_ptr_type(cpp_type) else ''
    f.write('  virtual void set_' + name + '(' + arg_type_const + get_type_as_ref_param(attr) + ' new_' + name + '_) {\n')
    # TODO: allow some values to be set even after initialization
    
    if has_single_superclass(class_def):
        f.write('    if (initialized) {\n')
        f.write('      throw RuntimeError("Value \'' + name + '\' of object with name " + name + " (class " + class_name + ") "\n')
        f.write('                         "cannot be set after model was initialized.");\n')
        f.write('    }\n')
        
    f.write('    ' + name + ' = new_' + name + '_;\n')
    f.write('  }\n') 

    # getter
    ret_type_const = 'const ' if is_cpp_ref_type(cpp_type) else ''
    
    f.write('  virtual ' + ret_type_const + get_type_as_ref_param(attr) + ' get_' + name + '() const {\n')
    f.write('    return ' + name + ';\n')
    f.write('  }\n') 
    return True
    
   
def write_method_signature(f, method):
    
    if KEY_RETURN_TYPE in method:
        f.write(yaml_type_to_cpp_type(method[KEY_RETURN_TYPE]) + ' ')
    else:
        f.write('void ')
    
    assert KEY_NAME in method
    f.write(rename_cpp_reserved_id(method[KEY_NAME]) + '(')
    
    if KEY_PARAMS in method:
        params = method[KEY_PARAMS]
        params_cnt = len(params)
        for i in range(params_cnt):
            p = params[i]
            assert KEY_NAME in p
            assert KEY_TYPE in p
            t = get_type_as_ref_param(p)
            const_spec = 'const ' if not is_yaml_ptr_type(p[KEY_TYPE]) and p[KEY_TYPE] != YAML_TYPE_PY_OBJECT else ''
            f.write(const_spec + t + ' ' + p[KEY_NAME])
            if KEY_DEFAULT in p:
                f.write(' = ' + get_default_or_unset_value(p))
            
            if i != params_cnt - 1:
                f.write(', ')
                
    f.write(')')
    
        
def write_method_declaration(f, method):
    f.write('  virtual ')
    write_method_signature(f, method)
    f.write(' = 0;\n')


def is_inherited(attr_or_method_def):
    if KEY_INHERITED in attr_or_method_def:
        return attr_or_method_def[KEY_INHERITED]
    else:
        return False
    

def attr_not_in_ctor(attr_or_method_def):
    if KEY_NOT_AS_CTOR_ARG in attr_or_method_def:
        return attr_or_method_def[KEY_NOT_AS_CTOR_ARG]
    else:
        return False
    

def has_single_superclass(class_def):
    if KEY_SUPERCLASS in class_def:
        return True
    else:
        return False


def has_superclass_other_than_base(class_def):
    # TODO: also deal with superclasses?
    if has_single_superclass(class_def):
        return \
            class_def[KEY_SUPERCLASS] != BASE_DATA_CLASS and \
            class_def[KEY_SUPERCLASS] != BASE_INTROSPECTION_CLASS
    else:
        return False
    
    
def write_gen_class(f, class_def, class_name):
    gen_class_name = GEN_CLASS_PREFIX + class_name
    f.write('class ' + gen_class_name)
    
    if has_single_superclass(class_def):
        f.write(': public ' + class_def[KEY_SUPERCLASS])
        
    f.write( ' {\n')
    f.write('public:\n')
    
    if has_superclass_other_than_base(class_def):
        write_ctor_for_superclass(f, class_def, class_name)
    
    if not has_single_superclass(class_def):
        # generate virtual destructor
        f.write('  virtual ~' + gen_class_name + '() {}\n')
        
    if has_single_superclass(class_def):
        f.write('  ' + RET_CTOR_POSTPROCESS + ' ' + CTOR_POSTPROCESS + '() ' + KEYWORD_OVERRIDE + ' {}\n')
        f.write('  ' + RET_TYPE_CHECK_SEMANTICS + ' ' + DECL_CHECK_SEMANTICS + ' ' + KEYWORD_OVERRIDE + ';\n')
        f.write('  void ' + DECL_SET_INITIALIZED + ' ' + KEYWORD_OVERRIDE + ';\n')
        # not virtual because overrides need different parameter type
        f.write('  bool __eq__(const ' + gen_class_name + '& other) const;\n')
        f.write('  ' +  RET_TYPE_SET_ALL_DEFAULT_OR_UNSET + ' ' + SET_ALL_DEFAULT_OR_UNSET_DECL + ' ' + KEYWORD_OVERRIDE + ';\n\n')

    override = KEYWORD_OVERRIDE if has_single_superclass(class_def) else ''
    f.write('  ' + RET_TYPE_TO_STR + ' ' + DECL_TO_STR_W_DEFAULT + ' ' + override + ';\n\n')
    
        
    f.write('  // --- attributes ---\n')
    items = class_def[KEY_ITEMS]
    for attr in items:
        if not is_inherited(attr):
            written = write_attr_with_get_set(f, class_def, attr)
            if written :
                f.write('\n')

    f.write('  // --- methods ---\n')
    methods = class_def[KEY_METHODS]
    for m in methods:
        if not is_inherited(m):
            write_method_declaration(f, m)
        
    f.write('}; // ' + GEN_CLASS_PREFIX + class_name + '\n\n')


def write_define_binding_decl(f, class_name, ret_type_void = False):
    if ret_type_void:
        f.write('void ')
    else:
        f.write('py::class_<' + class_name + '> ') 
    
    f.write('define_pybinding_' + class_name + '(py::module& m)')           


def remove_ptr_mark(t):
    assert len(t) > 1
    if t[-1] == '*':
        return t[0:-1]
    else:
        return t
    

def get_all_used_nonbase_types(class_def):
    types = set()
    for items in class_def[KEY_ITEMS]:
        assert KEY_TYPE in items
        t = items[KEY_TYPE]
        if not is_base_yaml_type(t):
            types.add( t )

    for method in class_def[KEY_METHODS]:
        if KEY_RETURN_TYPE in method:
            t = method[KEY_RETURN_TYPE]
            if not is_base_yaml_type(t):
                types.add( t )
            
        if KEY_PARAMS in method:
            for param in method[KEY_PARAMS]:
                assert KEY_TYPE in param
                t = param[KEY_TYPE]
                if not is_base_yaml_type(t):
                    types.add( t )

    return types


def get_all_used_compound_types_no_decorations(class_def):
    types = get_all_used_nonbase_types(class_def)
    cleaned_up_types = set()
    for t in types:
        # mus be set otherwise we get duplicates
        cleaned_up_types.add(remove_ptr_mark( get_first_inner_type( remove_ptr_mark(t) ) ))
    sorted_types_no_enums = [ t for t in cleaned_up_types if t not in g_enums ]
    sorted_types_no_enums.sort()
    return sorted_types_no_enums


def write_required_includes(f, class_def):
    types = get_all_used_nonbase_types(class_def)
    inner_types = set()
    for t in types:
        # mus be set otherwise we get duplicates
        inner_types.add( get_first_inner_type(t) )

    types_no_ptrs = [ t for t in inner_types if not is_yaml_ptr_type(t) ] 
    sorted_types_no_enums = [ t for t in types_no_ptrs if t not in g_enums ]
    sorted_types_no_enums.sort()
    
    for t in sorted_types_no_enums:
        f.write('#include "' + get_api_class_file_name_w_dir(t, EXT_H) + '"\n')


def write_forward_decls(f, class_def):
    # first we need to collect all types that we will need
    types = get_all_used_compound_types_no_decorations(class_def)
    
    for t in types:
        f.write('class ' + t + ';\n')
        
    if types:
        f.write('\n')

    
def generate_class_header(class_name, class_def):
    with open(get_gen_class_file_name_w_dir(class_name, EXT_H), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f)
        
        def_type = class_def[KEY_TYPE]
        
        guard = get_gen_header_guard_name(class_name);
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        f.write(INCLUDE_API_COMMON_H + '\n')
        
        if def_type == VALUE_CLASS:
            if KEY_SUPERCLASS in class_def:
                if class_def[KEY_SUPERCLASS] == BASE_DATA_CLASS:
                    f.write(INCLUDE_API_BASE_DATA_CLASS_H + '\n')
                elif class_def[KEY_SUPERCLASS] == BASE_INTROSPECTION_CLASS:
                    f.write(INCLUDE_API_BASE_INTROSPECTION_CLASS_H + '\n')
        
            if has_superclass_other_than_base(class_def):
                f.write('#include "' + get_api_class_file_name_w_dir(class_def[KEY_SUPERCLASS], EXT_H) + '"\n\n')
        
        write_required_includes(f, class_def)
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        write_forward_decls(f, class_def)
        
        if has_single_superclass(class_def): # not sure about this condition
            write_ctor_define(f, class_def, class_name)
        
        if def_type == VALUE_CLASS:
            write_gen_class(f, class_def, class_name)
        
            f.write('class ' + class_name + ';\n')
        elif def_type == VALUE_SUBMODULE:
            f.write('namespace ' + class_name + ' {\n\n')
            # submodules have only functions
            methods = class_def[KEY_METHODS]
            for m in methods:
                write_method_signature(f, m)
                f.write(";\n")
            f.write('\n} // namespace ' + class_name + '\n\n')
            
        write_define_binding_decl(f, class_name, def_type == VALUE_SUBMODULE)
        f.write(';\n')
        
        f.write(NAMESPACES_END + '\n\n')
        f.write('#endif // ' + guard + '\n')
    

def generate_class_template(class_name, class_def):
    with open(get_api_class_file_name_w_work_dir(class_name, EXT_H), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f)
        
        guard = get_api_header_guard_name(class_name);
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        f.write('#include "' + get_gen_class_file_name_w_dir(class_name, EXT_H) + '"\n')
        f.write(INCLUDE_API_COMMON_H + '\n')
        
        if has_superclass_other_than_base(class_def):
            f.write('#include "' + get_api_class_file_name_w_dir(class_def[KEY_SUPERCLASS], EXT_H) + '"\n\n')
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        
        f.write('class ' + class_name + ': public ' + GEN_CLASS_PREFIX + class_name + ' {\n')
        f.write('public:\n')
        f.write('  ' + get_underscored(class_name).upper() + CTOR_SUFFIX + '()\n')
        f.write('};\n\n')
        
        f.write(NAMESPACES_END + '\n\n')
        f.write('#endif // ' + guard + '\n')
        
        
def write_is_set_check(f, name):
    f.write('  if (!is_set(' + name + ')) {\n')
    f.write('    throw ValueError("Parameter \'' + name + '\' must be set.");\n')
    f.write('  }\n')
  

def write_check_semantics_implemetation(f, class_name, items):
    f.write(RET_TYPE_CHECK_SEMANTICS + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_CHECK_SEMANTICS + ' {\n') 
    for attr in items:
        if KEY_DEFAULT not in attr:
            write_is_set_check(f, attr[KEY_NAME])
    f.write('}\n\n')    
    
        
def write_to_str_implementation(f, class_name, items, based_on_base_superclass):
    f.write(RET_TYPE_TO_STR + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_TO_STR + ' {\n')

    # TODO: config etc from Model
    if class_name == "Model":
        f.write('  #if 0 // not generated correctly yet\n')

    f.write('  std::stringstream ss;\n')
    if based_on_base_superclass:
        f.write('  ss << get_object_name() << ": " <<\n')
    else:
        f.write('  ss << "' + class_name + '" << ": " <<\n')

    last_print_nl = False  
    
    num_attrs = len(items) 
    for i in range(num_attrs):
        name = items[i][KEY_NAME]
        type = items[i][KEY_TYPE]
        
        print_nl = False
        starting_nl = '"' if last_print_nl else '"\\n" << ind + "  " << "'
       
        if is_yaml_list_type(type):
            underlying_type = get_first_inner_type(type)
            if is_yaml_ptr_type(underlying_type):
                f.write('      ' + starting_nl + name + '=" << ' + VEC_PTR_TO_STR + '(' + name + ', ind + "  ")')
                print_nl = True
            else:
                f.write('      "' + name + '=" << ')
                f.write(VEC_NONPTR_TO_STR + '(' + name + ', ind + "  ")')
        elif not is_base_yaml_type(type) and type not in g_enums:
            f.write('      ' + starting_nl + name + '=" << "(" << ((' + name + ' != nullptr) ? ' + name + '->to_str(ind + "  ") : "null" ) << ")"')
            print_nl = True
        else:
            f.write('      "' + name + '=" << ')
            f.write(name)

        if i != num_attrs - 1:
            f.write(' << ", " <<')
            if print_nl:
                f.write(' "\\n" << ind + "  " <<\n')
            else:
                f.write('\n')
        
        last_print_nl = print_nl
        
    f.write(';\n')
    
    f.write('  return ss.str();\n')
    
    if class_name == "Model":
        f.write('  #else\n')
        f.write('  return "";\n')
        f.write('  #endif\n')
    
    f.write('}\n\n')                


def write_operator_equal_implemetation(f, class_name, items):
    # originally was operator== used, but this causes too cryptic compilation errors, 
    # so it was better to use python-style naming,
    # also the __eq__ is aonly defined for the 'Gen' classes, so defining operator ==
    # might have been more confusing 
    gen_class_name = GEN_CLASS_PREFIX + class_name
    f.write('bool ' + gen_class_name + '::__eq__(const ' + gen_class_name + '& other) const {\n')
    f.write('  return\n') 
    f.write('    name == other.name &&\n')
    
    num_attrs = len(items) 
    for i in range(num_attrs):
        name = items[i][KEY_NAME]
        t = items[i][KEY_TYPE]
        if is_yaml_ptr_type(t):
            
            f.write(
                '    (\n' 
                '      (' + name + ' != nullptr) ?\n'  
                '        ( (other.' + name + ' != nullptr) ?\n'
                '          (' + name + '->__eq__(*other.' + name + ')) : \n'
                '          false\n' 
                '        ) :\n'
                '        ( (other.' + name + ' != nullptr) ?\n'
                '          false :\n'
                '          true\n'
                '        )\n'
                '     ) '
            )
            
        elif is_yaml_list_type(t) and is_yaml_ptr_type(get_inner_list_type(t)):
            f.write('    vec_ptr_eq(' + name + ', other.' + name + ')')
        else:
            f.write('    ' + name + ' == other.' + name)
        if i != num_attrs - 1:
            f.write(' &&')
        else: 
            f.write(';')
        f.write('\n')
        
    f.write('}\n\n')    


def write_set_initialized_implemetation(f, class_name, items):
    # originally was operator== used, but this causes too cryptic compilation errors, 
    # so it was better to use python-style naming,
    # also the __eq__ is aonly defined for the 'Gen' classes, so defining operator ==
    # might have been more confusing 
    gen_class_name = GEN_CLASS_PREFIX + class_name
    f.write('void ' + gen_class_name + '::' + DECL_SET_INITIALIZED + ' {\n')
    num_attrs = len(items) 
    for i in range(num_attrs):
        name = items[i][KEY_NAME]
        t = items[i][KEY_TYPE]
        if is_yaml_ptr_type(t):
            f.write('  if (is_set(' + name + ')) {\n')
            f.write('    ' + name + '->set_initialized();\n')
            f.write('  }\n')
        elif is_yaml_list_type(t) and is_yaml_ptr_type(get_inner_list_type(t)):
            f.write('  vec_set_initialized(' + name + ');\n')
        
    f.write('  initialized = true;\n')    
    f.write('}\n\n')    
    
    
def is_overloaded(method, class_def):
    name = method[KEY_NAME] 
    count = 0
    for m in class_def[KEY_METHODS]:
        if m[KEY_NAME] == name:
            count += 1
             
    assert count >= 1
    return count >= 2 


def get_method_overload_cast(method):
    res = 'py::overload_cast<'
    
    if KEY_PARAMS in method:
        params = method[KEY_PARAMS]
        params_cnt = len(params)
        for i in range(params_cnt):
            p = params[i]
            assert KEY_TYPE in p
            t = get_type_as_ref_param(p)
            res += 'const ' + t
            if i + 1 != params_cnt:
                res += ', '
    
    res += '>'
    return res

def write_pybind11_method_bindings(f, class_name, method, class_def):
    assert KEY_NAME in method
    name = rename_cpp_reserved_id(method[KEY_NAME])
    
    full_method_name = '&' + class_name + '::' + name
    if is_overloaded(method, class_def):
        # overloaded method must be extra decorated with argument types for pybind11
        full_method_name = get_method_overload_cast(method) + '(' + full_method_name + ')' 
    
    f.write('      .def("' + name + '", ' + full_method_name)  
    
    if KEY_PARAMS in method:
        params = method[KEY_PARAMS]
        params_cnt = len(params)
        for i in range(params_cnt):
            f.write(', ')
            p = params[i]
            assert KEY_NAME in p
            f.write('py::arg("' + p[KEY_NAME] + '")')
            if KEY_DEFAULT in p:
                f.write(' = ' + get_default_or_unset_value(p))
                
    f.write(')\n')
    
    
def write_pybind11_bindings(f, class_name, class_def):
    def_type = class_def[KEY_TYPE]
    items = class_def[KEY_ITEMS]
    
    
    write_define_binding_decl(f, class_name, def_type != VALUE_CLASS)
    f.write(' {\n')
    
    superclass = ''
    if has_superclass_other_than_base(class_def):
        superclass = class_def[KEY_SUPERCLASS] + ', '

    ctor_has_args = True   
    if KEY_SUPERCLASS in class_def and class_def[KEY_SUPERCLASS] == BASE_INTROSPECTION_CLASS:
        ctor_has_args = False

    if def_type == VALUE_CLASS:
        f.write('  return py::class_<' + class_name + ', ' + superclass + \
                SHARED_PTR + '<' + class_name + '>>(m, "' + class_name + '")\n')
        f.write('      .def(\n')
        f.write('          py::init<\n')
    
        num_items = len(items)
        if ctor_has_args:
            if has_single_superclass(class_def):
                # init operands
                for i in range(num_items):
                    attr = items[i]
                    if attr_not_in_ctor(attr):
                        continue
                    const_spec = 'const ' if not is_yaml_ptr_type(attr[KEY_TYPE]) else ''
                    f.write('            ' + const_spec + get_type_as_ref_param(attr))
                    if i != num_items - 1:
                        f.write(',\n')
                if num_items != 0:
                    f.write('\n')
                
        f.write('          >()')
        
        if ctor_has_args:    
            # init argument names and default values
            if has_single_superclass(class_def):
                if num_items != 0:
                    f.write(',')
                f.write('\n')
            
                for i in range(num_items):
                    attr = items[i]
                    if attr_not_in_ctor(attr):
                        continue
                    name = attr[KEY_NAME]
                    f.write('          py::arg("' + name + '")')
                    if KEY_DEFAULT in attr:
                        f.write(' = ' + get_default_or_unset_value(attr))
                    if i != num_items - 1:
                        f.write(',\n')
        f.write('\n')          
        f.write('      )\n')            
        
        # common methods
        if has_single_superclass(class_def):
            f.write('      .def("check_semantics", &' + class_name + '::check_semantics)\n')
        
        f.write('      .def("__str__", &' + class_name + '::to_str, py::arg("ind") = std::string(""))\n')
        f.write('      .def("__repr__", &' + class_name + '::to_str, py::arg("ind") = std::string(""))\n')
    else:
        f.write('  m.def_submodule("' + class_name + '")\n')
        
    # declared methods
    for m in class_def[KEY_METHODS]:
        if not has_single_superclass(class_def) or not is_inherited(m):
            write_pybind11_method_bindings(f, class_name, m, class_def)

    if def_type == VALUE_CLASS:
        # dump needs to be always implemented
        f.write('      .def("dump", &' + class_name + '::dump)\n')
    
        # properties
        for i in range(num_items):
            if not has_single_superclass(class_def) or not is_inherited(items[i]):
                name = items[i][KEY_NAME]
                f.write('      .def_property("' + name + '", &' + class_name + '::get_' + name + ', &' + class_name + '::set_' + name + ')\n')
    f.write('    ;\n')
    f.write('}\n\n')
    
                
def write_used_classes_includes(f, class_def):
    types = get_all_used_compound_types_no_decorations(class_def)
    for t in types:
        f.write('#include "' + get_api_class_file_name_w_dir(t, EXT_H) + '"\n')


def write_set_all_default_or_unset(f, class_name, class_def):
    f.write(
        RET_TYPE_SET_ALL_DEFAULT_OR_UNSET + ' ' + GEN_CLASS_PREFIX + class_name + '::' + 
        SET_ALL_DEFAULT_OR_UNSET_DECL + ' {' + '\n'
    )
    
    # we do not need to call the derived methods because all members were inherited 
    
    f.write('  ' + CLASS_NAME_ATTR + ' = "' + class_name + '";\n')
    items = class_def[KEY_ITEMS]
    num_items = len(items)
    for i in range(num_items):
        assert KEY_NAME in items[i] 
        attr_name = items[i][KEY_NAME]
        f.write('  ' + attr_name + ' = ' + get_default_or_unset_value(items[i]) + ';\n')
        
    f.write('}\n\n')

            
def generate_class_implementation_and_bindings(class_name, class_def):
    with open(get_gen_class_file_name_w_dir(class_name, EXT_CPP), 'w') as f:
        def_type = class_def[KEY_TYPE]
        
        f.write(COPYRIGHT)
        write_generated_notice(f)
        
        f.write('#include <sstream>\n')
        f.write('#include "libs/pybind11/include/pybind11/stl.h"\n')

        # includes for our class
        f.write('#include "' + get_gen_class_file_name(class_name, EXT_H) + '"\n')
        
        if def_type == VALUE_CLASS:
            f.write('#include "' + get_api_class_file_name_w_dir(class_name, EXT_H) + '"\n')
        
        # we also need includes for every type that we used
        write_used_classes_includes(f, class_def)
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        if def_type == VALUE_CLASS:
            items = class_def[KEY_ITEMS]
            if has_single_superclass(class_def):
                write_check_semantics_implemetation(f, class_name, items)
                write_operator_equal_implemetation(f, class_name, items)
                write_set_initialized_implemetation(f, class_name, items)
                write_set_all_default_or_unset(f, class_name, class_def)
    
            write_to_str_implementation(f, class_name, items, has_single_superclass(class_def))
        
        write_pybind11_bindings(f, class_name, class_def)
        
        f.write(NAMESPACES_END + '\n\n')
        

# class Model provides the same methods as Subsystem and InstantiationData, 
# this function copies the definition for pybind11 API generation
def inherit_from_superclasses(data_classes, class_name, class_def):
    
    res = class_def.copy() 
    
    superclass_names = []
    if has_superclass_other_than_base(class_def):
        superclass_names.append(class_def[KEY_SUPERCLASS])
        
    if KEY_SUPERCLASSES in class_def:
        superclass_names += class_def[KEY_SUPERCLASSES]
        
    for sc_name in superclass_names:
        assert sc_name in data_classes
        superclass_def = data_classes[sc_name]
        assert not has_superclass_other_than_base(superclass_def), "Only one level of inheritance"
        
        # we are not checking any duplicates
        if KEY_METHODS in superclass_def:
            for method in superclass_def[KEY_METHODS]:
                method_copy = copy(method)
                method_copy[KEY_INHERITED] = True
                res[KEY_METHODS].append(method_copy) 
        
        class_attr_names = set()
        if KEY_ITEMS in class_def:
            for attr in class_def[KEY_ITEMS]:
                class_attr_names.add(attr[KEY_NAME])

        if KEY_ITEMS in superclass_def:
            for attr in superclass_def[KEY_ITEMS]:
                
                attr_copy = copy(attr)
                attr_copy[KEY_INHERITED] = True
                res[KEY_ITEMS].append(attr_copy) 
                
                # do not inherit attribute if the attribute with the same name already exists
                if attr[KEY_NAME] in class_attr_names:
                    attr_copy[KEY_NOT_AS_CTOR_ARG] = True               
        
    return res


def generate_class_files(data_classes, class_name, class_def):
    # we need items and methods to be present 
    if KEY_ITEMS not in class_def:
        class_def[KEY_ITEMS] = []
        
    if KEY_METHODS not in class_def:
        class_def[KEY_METHODS] = []

    class_def_w_inheritances = inherit_from_superclasses(data_classes, class_name, class_def)
    
    generate_class_header(class_name, class_def_w_inheritances)
    generate_class_implementation_and_bindings(class_name, class_def_w_inheritances)
    
    generate_class_template(class_name, class_def_w_inheritances)
    
    
def write_constant_def(f, constant_def):
    assert KEY_NAME in constant_def
    assert KEY_TYPE in constant_def
    assert KEY_VALUE in constant_def
    name = constant_def[KEY_NAME]
    t = constant_def[KEY_TYPE]
    value = constant_def[KEY_VALUE]
    
    # complex types are defined manually in globals.h
    if is_base_yaml_type(t):
        q = '"' if t == YAML_TYPE_STR else ''
        f.write('const ' + yaml_type_to_cpp_type(t) + ' ' + name + ' = ' + q + str(value) + q + ';\n')
        
        
def write_enum_def(f, enum_def):
    assert KEY_NAME in enum_def
    assert KEY_VALUES in enum_def, "Enum must have at least one value"
    name = enum_def[KEY_NAME]
    values = enum_def[KEY_VALUES]
    
    f.write('\nenum class ' + name + ' {\n')
    num = len(values)
    for i in range(num):
        v = values[i]
        assert KEY_NAME in v
        assert KEY_VALUE in v
        assert type(v[KEY_VALUE]) == int  
        f.write('  ' + v[KEY_NAME] + ' = ' + str(v[KEY_VALUE]))
        if i + 1 != num:
            f.write(',')
        f.write('\n')    
    f.write('};\n\n')

    
    f.write('\nstatic inline  std::ostream& operator << (std::ostream& out, const ' + name + ' v) {\n')
    f.write('  switch (v) {\n');
    for i in range(num):
        v = values[i]
        # dumping in Python style '.' that will be most probably more common as API  
        f.write('    case ' + name + '::' + v[KEY_NAME] + ': out << "' + name + '.' + v[KEY_NAME] + ' (' + str(v[KEY_VALUE]) + ')"; break;\n')  
        
    f.write('  }\n')
    f.write('  return out;\n')
    f.write('};\n\n')
        
            
def generate_constants_header(constants_items, enums_items):
    with open(os.path.join(TARGET_DIRECTORY, GEN_CONSTANTS_H), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f)
        
        guard = 'API_GEN_CONSTANTS';
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        
        f.write('#include <string>\n')
        f.write('#include "api/globals.h"\n')
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        for constant_def in constants_items:
            write_constant_def(f, constant_def)
            
        for enum_def in enums_items:
            write_enum_def(f, enum_def)

        f.write('\n' + DECL_DEFINE_PYBINDIND_CONSTANTS + ';\n')
        f.write(DECL_DEFINE_PYBINDIND_ENUMS + ';\n\n')
        
        f.write(NAMESPACES_END + '\n\n')
        
        f.write('#endif // ' + guard + '\n\n')


def write_constant_binding(f, constant_def):
    assert KEY_NAME in constant_def
    assert KEY_TYPE in constant_def
    assert KEY_VALUE in constant_def
    name = constant_def[KEY_NAME]
    t = constant_def[KEY_TYPE]
    value = constant_def[KEY_VALUE]
    
    #q = '"' if t == YAML_TYPE_STR else ''
    
    pybind_type = yaml_type_to_pybind_type(t)
    if pybind_type == PYBIND_TYPE_OBJECT:
        f.write('  m.attr("' + name + '") = py::' + PYBIND_TYPE_OBJECT + '(' + PY_CAST + '(' + name + '));\n')
    else:
        f.write('  m.attr("' + name + '") = py::' + pybind_type + '(' + name + ');\n')
    

def write_enum_binding(f, enum_def):
    assert KEY_NAME in enum_def
    assert KEY_VALUES in enum_def, "Enum must have at least one value"
    name = enum_def[KEY_NAME]
    values = enum_def[KEY_VALUES]
    
    f.write('  py::enum_<' + name + '>(m, "' + name + '", py::arithmetic())\n')
    num = len(values)
    for i in range(num):
        v = values[i]
        assert KEY_NAME in v
        f.write('    .value("' + v[KEY_NAME] + '", ' + name + '::' + v[KEY_NAME] + ')\n')
        
    f.write('    .export_values();\n')
        
        
def generate_constants_implementation(constants_items, enums_items):
    with open(os.path.join(TARGET_DIRECTORY, GEN_CONSTANTS_CPP), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f)
        f.write(INCLUDE_API_COMMON_H + '\n')
        
        # includes for used classes
        used_classes = set()
        for constant_def in constants_items:
            t = constant_def[KEY_TYPE]
            pybind_type = yaml_type_to_pybind_type(t)
            if pybind_type == PYBIND_TYPE_OBJECT:
                used_classes.add(t)
        for cls in used_classes:
            f.write('#include "' + get_api_class_file_name_w_dir(cls, EXT_H) + '"\n')
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        f.write(DECL_DEFINE_PYBINDIND_CONSTANTS + ' {\n')
        
        for constant_def in constants_items:
            write_constant_binding(f, constant_def)

        f.write('}\n\n')

        f.write(DECL_DEFINE_PYBINDIND_ENUMS + ' {\n')
        
        for enum_def in enums_items:
            write_enum_binding(f, enum_def)
        
        f.write('}\n\n')
        
        f.write(NAMESPACES_END + '\n\n')    


def generate_constants_and_enums(constants_items, enums_items):
    generate_constants_header(constants_items, enums_items)
    generate_constants_implementation(constants_items, enums_items)
    

def set_global_enums(data_classes):
    global g_enums
    res = set()
    enum_defs = data_classes[KEY_ENUMS] if KEY_ENUMS in data_classes else []
    for enum in enum_defs:
        assert KEY_NAME in enum
        res.add(enum[KEY_NAME])
    g_enums = res


def generate_data_classes(data_classes):
    generate_constants_and_enums(
        data_classes[KEY_CONSTANTS] if KEY_CONSTANTS in data_classes else [],
        data_classes[KEY_ENUMS] if KEY_ENUMS in data_classes else [])

    set_global_enums(data_classes)

    for key, value in data_classes.items():
        if key != KEY_CONSTANTS and key != KEY_ENUMS:
            # set default definition type
            if not KEY_TYPE in value:
                value[KEY_TYPE] = VALUE_CLASS
            
            if VERBOSE:
                if value[KEY_TYPE] == VALUE_CLASS:
                    print("Generating class " + key)
                elif value[KEY_TYPE] == VALUE_SUBMODULE:
                    print("Generating submodule " + key)
                else:
                    assert False, "Unknown definition type"     
            
            generate_class_files(data_classes, key, value)
    
    
def collect_all_names(data_classes):
    all_class_names = set()    
    all_item_param_names = set()
    all_enum_value_names = set()
    all_const_value_names = set()
    
    for key, value in data_classes.items():
        
        if key != KEY_CONSTANTS and key != KEY_ENUMS:
            all_class_names.add(key)
            # items 
            if KEY_ITEMS in value:
                for item in value[KEY_ITEMS]:
                    all_item_param_names.add(item[KEY_NAME])
                    
            # methods        
            if KEY_METHODS in value:
                for method in value[KEY_METHODS]:
                    all_item_param_names.add(method[KEY_NAME])

                    if KEY_PARAMS in method:
                        for param in method[KEY_PARAMS]:
                            all_item_param_names.add(param[KEY_NAME])
        elif key == KEY_ENUMS:
            # enum names are in g_enums
            for enum in value:
                if KEY_VALUES in enum:
                    for enum_item in enum[KEY_VALUES]:
                        all_enum_value_names.add(enum_item[KEY_NAME])
        elif key == KEY_CONSTANTS:
            for const in value:
                all_const_value_names.add(const[KEY_NAME])
        else:
            print("Error: unexpected top level key " + key)
            sys.exit(1)
    
    all_class_names = list(all_class_names)
    all_class_names.sort(key=str.casefold)
    
    all_item_param_names_list = list(all_item_param_names)
    all_item_param_names_list.sort(key=str.casefold)
    
    all_enum_value_names_list = list(all_enum_value_names)
    all_enum_value_names_list.sort(key=str.casefold)
    
    all_const_value_names_list = list(all_const_value_names)
    all_const_value_names_list.sort(key=str.casefold)
    
    return all_class_names, all_item_param_names_list, all_enum_value_names_list, all_const_value_names_list

  
def write_name_def(f, name, extra_prefix=''):
    upper_name = get_underscored(name).upper()
    f.write('const char* const ' + NAME_PREFIX + extra_prefix + upper_name + ' = "' + name + '";\n')


def write_name_def_verbatim(f, name, extra_prefix=''):
    f.write('const char* const ' + NAME_PREFIX + extra_prefix + name + ' = "' + name + '";\n')
      
      
# this function generates definitions for converter so that we can use constant strings 
def generate_names_header(data_classes):

    all_class_names, all_item_param_names_list, all_enum_value_names_list, all_const_value_names_list = collect_all_names(data_classes)
    
    with open(os.path.join(TARGET_DIRECTORY, GEN_NAMES_H), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f)
        
        guard = 'API_GEN_NAMES';
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')

        for name in all_class_names:
            write_name_def(f, name, CLASS_PREFIX)
        f.write('\n')

        for name in all_item_param_names_list:
            write_name_def(f, name)
        f.write('\n')

        all_enums_list = list(g_enums)
        all_enums_list.sort(key=str.casefold)
        for name in all_enums_list:
            write_name_def(f, name, ENUM_PREFIX)
        f.write('\n')
        
        for name in all_enum_value_names_list:
            write_name_def_verbatim(f, name, ENUM_VALUE_PREFIX)
        f.write('\n')
        
        for name in all_const_value_names_list:
            write_name_def_verbatim(f, name, CONSTANT_VALUE_PREFIX)
        f.write('\n')
        
        f.write(NAMESPACES_END + '\n\n')
        f.write('#endif // ' + guard + '\n\n')      
    
    
def generate_pyi_class(f, name, class_def):
    f.write('class ' + name + '():\n')
    f.write('    def __init__(\n')
    param_ind = '            ' 
    f.write(param_ind + 'self,\n')
    
    # ctor
    if KEY_ITEMS in class_def:
        generated_params = set()
        num_items = len(class_def[KEY_ITEMS])
        for i in range(0, num_items):
            item = class_def[KEY_ITEMS][i]
            name = item[KEY_NAME]
            if name in generated_params:
                continue
            generated_params.add(name)
            
            f.write(param_ind + name)
            
            # self-referencing classes must use string type reference
            # https://www.python.org/dev/peps/pep-0484/#forward-references
            t = yaml_type_to_py_type(item[KEY_TYPE])
            q = '\'' if t == name else ''
            f.write(' : ' + q + t + q)
            
            if KEY_DEFAULT in item:
                f.write(' = ' + get_default_or_unset_value_py(item))
            if i + 1 != num_items:
                f.write(',')
            f.write('\n')
    f.write('        ):\n')

    # class members
    if KEY_ITEMS in class_def and class_def[KEY_ITEMS]:
        num_items = len(class_def[KEY_ITEMS])
        for i in range(0, num_items):
            member_name = class_def[KEY_ITEMS][i][KEY_NAME]
            f.write('        self.' + member_name + ' = ' + member_name + '\n')
    else:
        f.write('        pass')
    f.write('\n\n')        
        
        
    if KEY_METHODS in class_def:
        # for now, we just simply print the first variant oif there are multiple methods with the same name
        printed_methods = set()
        
        for method in class_def[KEY_METHODS]:
            method_name = method[KEY_NAME]
            if method_name in printed_methods:
                continue
            printed_methods.add(method_name)
            
            f.write('    def ' + method[KEY_NAME] + '(\n')
            f.write(param_ind + 'self,\n')
            
            if KEY_PARAMS in method:
                num_params = len(method[KEY_PARAMS])
                for i in range(0, num_params):
                    param = method[KEY_PARAMS][i]
                    f.write(param_ind + param[KEY_NAME])
                    t = yaml_type_to_py_type(param[KEY_TYPE])
                    q = '\'' if t == name else ''
                    f.write(' : ' + q + t + q)
                    if KEY_DEFAULT in param:
                        f.write(' = ' + get_default_or_unset_value_py(param))
                    if i + 1 != num_params:
                        f.write(',')       
                    f.write('\n')   
            f.write('        )')
            
            if KEY_RETURN_TYPE in method:
                f.write(' -> \'' + yaml_type_to_py_type(method[KEY_RETURN_TYPE]) + '\'')
            else:
                f.write(' -> ' + PY_NONE)
            f.write(':\n')
            f.write('        pass\n\n')
            
                
def generate_pyi_file(data_classes):
    with open(os.path.join(TARGET_DIRECTORY, MCELL_PYI), 'w') as f:
        
        f.write('from typing import List, Dict, Callable, Any\n')
        f.write('from enum import Enum\n\n')
        
        f.write('# "forward" declarations to make the type hints valid\n')
        for key, value in sorted(data_classes.items()):
            if key != KEY_CONSTANTS and key != KEY_ENUMS:
                f.write('class ' + key + '():\n    pass\n')
        f.write('\n')
            
        species_def = ''
        
        # Vec3
        f.write(
            'class Vec3():\n'
            '    def __init__(self, x : float = 0, y : float = 0, z : float = 0):\n'
            '        self.x = x\n'
            '        self.y = y\n'
            '        self.z = z\n'
            '\n'
        )
        
        # generate enums first, then constants
        enums = data_classes[KEY_ENUMS]
        for enum in enums:
            f.write('class ' + enum[KEY_NAME] + '(Enum):\n')
            for enum_item in enum[KEY_VALUES]:
                f.write('    ' + enum_item[KEY_NAME] + ' = ' + str(enum_item[KEY_VALUE]) + '\n')
            f.write('\n')
        f.write('\n\n')
        
        constants = data_classes[KEY_CONSTANTS]
        for const in constants:
            if const[KEY_TYPE] == YAML_TYPE_SPECIES:
                # Species constants are a bit special
                ctor_param = get_underscored(const[KEY_VALUE]).upper()
                species_def += const[KEY_NAME] + ' = ' + YAML_TYPE_SPECIES + '(\'' + ctor_param + '\')\n'
            elif const[KEY_TYPE] == YAML_TYPE_STR:
                f.write(const[KEY_NAME] + ' = \'' + str(const[KEY_VALUE]) + '\'\n')
            else:
                f.write(const[KEY_NAME] + ' = ' + str(const[KEY_VALUE]) + '\n')
        f.write('\n\n')
            
        for key, value in sorted(data_classes.items()):
            if key != KEY_CONSTANTS and key != KEY_ENUMS:
                generate_pyi_class(f, key, value)

        # we need to define the Species contants after they were defined
        f.write(species_def)
        
        
def indent_and_fix_rst_chars(text, indent_and_fix_rst_chars):
    return text.replace('\n', '\n' + indent_and_fix_rst_chars).replace('*', '\\*')
            

def generate_class_documentation(f, class_name, class_def):
    f.write(class_name + '\n' + '='*len(class_name) + '\n\n')
    
    if KEY_DOC in class_def:
        f.write(class_def[KEY_DOC].strip() + '\n\n')
        
    if KEY_ITEMS in class_def and class_def[KEY_ITEMS]:
        f.write('Attributes:\n' + '*'*len('Attributes:') + '\n')
        num_items = len(class_def[KEY_ITEMS])
        for item in class_def[KEY_ITEMS]:
            t = yaml_type_to_py_type(item[KEY_TYPE])
            f.write('* | **' + item[KEY_NAME] + '**: ' + t)
            
            if KEY_DEFAULT in item:
                f.write(' = ' + get_default_or_unset_value_py(item))
            f.write('\n')
            
            if KEY_DOC in item:
                f.write('  | ' + indent_and_fix_rst_chars(item[KEY_DOC].strip(), '  | ') + '\n')
            f.write('\n')
        
    if KEY_METHODS in class_def and class_def[KEY_METHODS]:
        f.write('\nMethods:\n' + '*'*len('nMethods:') + '\n')
        for method in class_def[KEY_METHODS]:
            method_name = method[KEY_NAME]
            
            f.write('* | **' + method[KEY_NAME] + '**\n\n')

            if KEY_PARAMS in method:
                num_params = len(method[KEY_PARAMS])
                for param in method[KEY_PARAMS]:
                    t = yaml_type_to_py_type(param[KEY_TYPE])
                    f.write('   * | ' + param[KEY_NAME] + ': ' + t)
                    if KEY_DEFAULT in param:
                        f.write(' = ' + get_default_or_unset_value_py(param))
                    if KEY_DOC in param:
                        f.write('\n     | ' + indent_and_fix_rst_chars(param[KEY_DOC].strip(), '     | ') + '\n\n')
                    else:
                        f.write('\n')

            if KEY_RETURN_TYPE in method:
                f.write('   * | return type: ' + yaml_type_to_py_type(method[KEY_RETURN_TYPE]) + '\n\n')
            
            if KEY_DOC in method:
                f.write('\n  | ' + indent_and_fix_rst_chars(method[KEY_DOC].strip(), '  | ') + '\n\n')
                


            f.write('\n')
        f.write('\n')
    
        
def generate_documentation(data_classes):
    with open(os.path.join(DOC_DIRECTORY, CLASSES_RST), 'w') as f:

        f.write(
            '**************************************\n'
            'MCell 4 Python API Enums and Constants\n'
            '**************************************\n\n'
        )

        # generate enums first, then constants
        enums = data_classes[KEY_ENUMS]
        for enum in enums:
            enum_name = enum[KEY_NAME]
            f.write(enum_name + '\n' + '='*len(enum_name) + '\n\n')
            
            if KEY_DOC in enum:
                f.write('\n  | ' + indent_and_fix_rst_chars(enum[KEY_DOC].strip(), '  | ') + '\n\n')            
                
            for value in enum[KEY_VALUES]:
                f.write('* | **' + value[KEY_NAME] + '** = ' + str(value[KEY_VALUE]) + '\n')
                if KEY_DOC in value:
                    f.write('  | ' + indent_and_fix_rst_chars(value[KEY_DOC].strip(), '  | ') + '\n\n')                
            f.write('\n')
        f.write('\n\n')
        
        c = 'Constants'
        f.write(c + '\n' + '='*len(c) + '\n\n')
        constants = data_classes[KEY_CONSTANTS]
        for const in constants:
            const_name = const[KEY_NAME]
            
            f.write('* | **' + const_name + '**: ' + yaml_type_to_py_type(const[KEY_TYPE]) + \
                    ' = ' + str(const[KEY_VALUE]) +'\n')
            if KEY_DOC in const:
                f.write('  | ' + indent_and_fix_rst_chars(const[KEY_DOC].strip(), '  | ') + '\n\n')              
        f.write('\n\n')
        
        f.write(
            '**************************\n'
            'MCell 4 Python API Classes\n'
            '**************************\n\n'
        )
                    
        for key, value in sorted(data_classes.items()):
            if key != KEY_CONSTANTS and key != KEY_ENUMS:
                generate_class_documentation(f, key, value)


        
def load_and_generate_data_classes():
    # create work directory
    work_dir = os.path.join('..', 'work')
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    
    data_classes = {}
    
    for input_file in ALL_INPUT_FILES:
        with open(input_file) as file:
            # The FullLoader parameter handles the conversion from YAML
            # scalar values to Python the dictionary format
            loaded_data = yaml.load(file, Loader=yaml.FullLoader)
            if loaded_data:
                data_classes.update(loaded_data)
            else:
                print("File " + input_file + " is empty (or could not be loaded as yaml), ignoring it")
            if VERBOSE:
                print("Loaded " + input_file)
    
    if VERBOSE:    
        print(data_classes)
    assert type(data_classes) == dict
    generate_data_classes(data_classes)
    generate_names_header(data_classes)
    generate_pyi_file(data_classes)
    generate_documentation(data_classes)

if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == '-v':
        VERBOSE = True
        
    load_and_generate_data_classes()
    
    