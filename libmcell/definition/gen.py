#!/usr/bin/env python3
"""
Copyright (C) 2020 by
The Salk Institute for Biological Studies 

Use of this source code is governed by an MIT-style
license that can be found in the LICENSE file or at
https://opensource.org/licenses/MIT.
"""

# TODO: change 'items' to 'attributes'?
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
from pyexpat import model

VERBOSE = False # may be overridden by argument -v 

from constants import *
import doc

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
    return GEN_PREFIX + get_underscored(class_name) + '.' + extension

def get_gen_class_file_name_w_dir(class_name, extension):
    # using '/' intentionally to generate the same output on Windows and Linux/Mac 
    return TARGET_DIRECTORY + '/' + get_gen_class_file_name(class_name, extension)

def get_api_class_file_name(class_name, extension):
    return get_underscored(class_name) + '.' + extension

def get_api_class_file_name_w_dir(class_name, extension):
    if class_name == "Config":
        # config.h is be called api_config.h due to include collisions with MSVC 
        class_name = 'ApiConfig'

    return API_DIRECTORY + '/' + get_api_class_file_name(class_name, extension)
   
def get_api_class_file_name_w_work_dir(class_name, extension):
    return WORK_DIRECTORY + '/' + get_api_class_file_name(class_name, extension)

def get_as_shared_ptr(class_name):
    return SHARED_PTR + '<' + class_name + '>'

def get_copy_function_name(class_name):
    return COPY_NAME + '_' + get_underscored(class_name)

def get_deepcopy_function_name(class_name):
    return DEEPCOPY_NAME + '_' + get_underscored(class_name)

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

# returns true also for enums
def is_base_yaml_type(t):
    return \
        t == YAML_TYPE_FLOAT or t == YAML_TYPE_STR or t == YAML_TYPE_INT or t == YAML_TYPE_UINT64 or t == YAML_TYPE_UINT32 or \
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


def yaml_type_to_cpp_type(t, w_namespace=False):
    assert len(t) >= 1
    if t == YAML_TYPE_FLOAT:
        return CPP_TYPE_DOUBLE
    elif t == YAML_TYPE_STR:
        return CPP_TYPE_STR
    elif t == YAML_TYPE_INT:
        return CPP_TYPE_INT
    elif t == YAML_TYPE_UINT64:
        return CPP_TYPE_UINT64
    elif t == YAML_TYPE_UINT32:
        return CPP_TYPE_UINT32
    elif t == YAML_TYPE_BOOL:
        return CPP_TYPE_BOOL
    elif t == YAML_TYPE_VEC2:
        ns = '' if not w_namespace else 'MCell::' 
        return ns + CPP_TYPE_VEC2
    elif t == YAML_TYPE_VEC3:
        ns = '' if not w_namespace else 'MCell::' 
        return ns + CPP_TYPE_VEC3
    elif t == YAML_TYPE_IVEC3:
        ns = '' if not w_namespace else 'MCell::' 
        return ns + CPP_TYPE_IVEC3
    elif is_yaml_list_type(t):
        assert len(t) > 7
        inner_type = yaml_type_to_cpp_type(get_inner_list_type(t), w_namespace)
        return CPP_VECTOR_TYPE + '<' + inner_type + '>'
    elif is_yaml_dict_type(t):
        key_type = yaml_type_to_cpp_type(get_inner_dict_key_type(t), w_namespace)
        value_type = yaml_type_to_cpp_type(get_inner_dict_value_type(t), w_namespace)
        return CPP_MAP_TYPE + '<' + key_type + ', ' + value_type + '>'
    else:
        ns = '' if not w_namespace else 'MCell::API::' 
        if is_yaml_ptr_type(t):
            return SHARED_PTR + '<' + ns + t[0:-1] + '>' 
        else:
            return ns + t # standard ttype
    
    
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
    elif t == YAML_TYPE_UINT64:
        return CPP_TYPE_INT + '_'
    elif t == YAML_TYPE_UINT32:
        return CPP_TYPE_INT + '_'
    elif t == YAML_TYPE_SPECIES:
        return PYBIND_TYPE_OBJECT
    else:
        assert False, "Unsupported constant type " + t 
    
    
def yaml_type_to_py_type(t):
    assert len(t) >= 1
    if t == YAML_TYPE_UINT64 or t == YAML_TYPE_UINT32:
        return YAML_TYPE_INT # not sure what should be the name
    elif is_yaml_function_type(t):
        return "Callable, # " + t
    elif t == YAML_TYPE_PY_OBJECT:
        return "Any, # " + t
    elif is_yaml_list_type(t) and get_inner_list_type(t) == YAML_TYPE_UINT64:
        return t.replace(YAML_TYPE_UINT64, YAML_TYPE_INT)
    elif is_yaml_list_type(t) and get_inner_list_type(t) == YAML_TYPE_UINT32:
        return t.replace(YAML_TYPE_UINT32, YAML_TYPE_INT)
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
    
def is_cpp_vector_type(cpp_type):
    return'std::vector' in cpp_type
    
        
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
    elif t == YAML_TYPE_UINT64:
        return UNSET_VALUE_UINT64
    elif t == YAML_TYPE_UINT32:
        return UNSET_VALUE_UINT32
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


def has_item_w_name(items, item_name):
    for item in items:
        if item_name == item[KEY_NAME]:
            return True
    return False


def is_container_class(name):
    return \
        name == CLASS_NAME_MODEL or \
        name == CLASS_NAME_SUBSYSTEM or \
        name == CLASS_NAME_INSTANTIATION or \
        name == CLASS_NAME_OBSERVABLES

def is_container_class_no_model(name):
    return is_container_class(name) and not name == CLASS_NAME_MODEL

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
            # and we need to initialize Complex 
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
    

def needs_default_ctor(class_def, only_inherited):
    
    if KEY_ITEMS not in class_def:
        return False
    
    items = class_def[KEY_ITEMS]
    if only_inherited:
        items = [ attr for attr in items if is_inherited(attr) ]
        
    if items:
        all_attrs_initialized = True
        for item in items:
            if not KEY_DEFAULT in item:
                 all_attrs_initialized = False
                 
        return not all_attrs_initialized
    else:
        return False
        

def write_ctor_define(f, class_def, class_name):
    with_args = True
    if KEY_SUPERCLASS in class_def and class_def[KEY_SUPERCLASS] == BASE_INTROSPECTION_CLASS:
        with_args = False
        
    suffix = CTOR_SUFFIX if with_args else CTOR_NOARGS_SUFFIX
    f.write('#define ' + get_underscored(class_name).upper() + suffix + '() \\\n')

    write_ctor_decl(f, class_def, class_name, append_backslash=True, indent_and_fix_rst_chars='    ', only_inherited=False, with_args=with_args)

    f.write('{ \\\n')

    # initialization code
    if not is_container_class_no_model(class_name):
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
            
    if not is_container_class_no_model(class_name):
        f.write('      ' + CTOR_POSTPROCESS + '(); \\\n')    
        f.write('      ' + CHECK_SEMANTICS + '(); \\\n')
    f.write('    } \\\n')
    
    # also generate empty ctor if needed
    f.write('    ' + class_name + '(' + DEFAULT_CTOR_ARG_TYPE + ')')
    
    if has_single_superclass(class_def):
        f.write(' : \\\n' 
                '      ' + GEN_CLASS_PREFIX + class_name + '(' + DEFAULT_CTOR_ARG_TYPE + '()) ')
        
    f.write('{ \\\n')
    
    if has_single_superclass(class_def):
        f.write('      ' + SET_ALL_DEFAULT_OR_UNSET_DECL + '; \\\n')
        f.write('      ' + SET_ALL_CUSTOM_TO_DEFAULT_DECL + '; \\\n')
    
    f.write('    }\n\n');
    
    
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
    
    if has_single_superclass(class_def):
        f.write('    ' + CACHED_DATA_ARE_UPTODATE_ATTR + ' = false;\n');
    f.write('    ' + name + ' = new_' + name + '_;\n')
    f.write('  }\n') 

    # getter
    ret_type_ref = ''
    ret_type_const = 'const ' if is_cpp_ref_type(cpp_type) else ''
    method_const = 'const '
    
    # vectors are always returned as a non-const reference
    if is_cpp_vector_type(cpp_type):
        ret_type_ref = '&'
        ret_type_const = ''
        method_const = ''
    
    ret_type = get_type_as_ref_param(attr)
    f.write('  virtual ' + ret_type_const + ret_type + ret_type_ref + ' get_' + name + '() ' + method_const +'{\n')
    """
    # original approach to protect already initialized classes but this behavior woudl probably be more confusing
    if has_single_superclass(class_def) and is_cpp_vector_type(cpp_type):
        f.write('    if (initialized) {\n')
        f.write('      // any changes after initialization are ignored\n')
        f.write('      static ' + ret_type + ' copy = ' + name + ';\n')
        f.write('      return copy;\n')
        f.write('    }\n')
    """
    if has_single_superclass(class_def):
        f.write('    ' + CACHED_DATA_ARE_UPTODATE_ATTR + ' = false; // arrays and other data can be modified through getters\n');
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
    
    if KEY_IS_CONST in method and method[KEY_IS_CONST]:
        f.write(' const')
    
        
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
    
    
def write_gen_class(f, class_def, class_name, decls):
    gen_class_name = GEN_CLASS_PREFIX + class_name
    f.write('class ' + gen_class_name)
    
    if has_single_superclass(class_def):
        f.write(': public ' + class_def[KEY_SUPERCLASS])
    elif KEY_SUPERCLASSES in class_def:
        scs = class_def[KEY_SUPERCLASSES]
        assert scs
        f.write(': ')
        for i in range(len(scs)):
            f.write('public ' + scs[i])
            if i + 1 != len(scs):
                f.write(', ')
    elif is_container_class_no_model(class_name):
        f.write(': public ' + BASE_EXPORT_CLASS)
        
    f.write( ' {\n')
    f.write('public:\n')
    
    initialization_ctor_generated = False
    if has_superclass_other_than_base(class_def):
        write_ctor_for_superclass(f, class_def, class_name)
        initialization_ctor_generated = True
    
    # we need a 'default' ctor (with custom arg type) all the time,
    # however, there may not be a default ctor with no argument so we must define it 
    # as well
    if needs_default_ctor(class_def, only_inherited=True) or not initialization_ctor_generated:
        
        f.write('  ' + GEN_CLASS_PREFIX + class_name + '() ')
        
        if has_superclass_other_than_base(class_def):
            superclass_name = class_def[KEY_SUPERCLASS]
            f.write(': ' + superclass_name + '(' + DEFAULT_CTOR_ARG_TYPE + '()) {\n')
        else:
            f.write('{\n')
        f.write('  }\n')
    
    f.write('  ' + GEN_CLASS_PREFIX + class_name + '(' + DEFAULT_CTOR_ARG_TYPE + ') {\n')
    f.write('  }\n')


    if not has_single_superclass(class_def):
        # generate virtual destructor
        f.write('  virtual ~' + gen_class_name + '() {}\n')
        
    if has_single_superclass(class_def):
        f.write('  ' + RET_CTOR_POSTPROCESS + ' ' + CTOR_POSTPROCESS + '() ' + KEYWORD_OVERRIDE + ' {}\n')
        f.write('  ' + RET_TYPE_CHECK_SEMANTICS + ' ' + DECL_CHECK_SEMANTICS + ' ' + KEYWORD_OVERRIDE + ';\n')
        f.write('  void ' + DECL_SET_INITIALIZED + ' ' + KEYWORD_OVERRIDE + ';\n')
        f.write('  ' +  RET_TYPE_SET_ALL_DEFAULT_OR_UNSET + ' ' + SET_ALL_DEFAULT_OR_UNSET_DECL + ' ' + KEYWORD_OVERRIDE + ';\n\n')

    f.write('  ' + get_as_shared_ptr(class_name) + ' ' + get_copy_function_name(class_name) + '() const;\n')
    f.write('  ' + get_as_shared_ptr(class_name) + ' ' + get_deepcopy_function_name(class_name) + '(py::dict = py::dict()) const;\n')
    f.write('  virtual bool __eq__(const ' + class_name + '& other) const;\n')
    f.write('  virtual bool eq_nonarray_attributes(const ' + class_name + '& other, const bool ignore_name = false) const;\n')
    f.write('  bool operator == (const ' + class_name + '& other) const { return __eq__(other);}\n')
    f.write('  bool operator != (const ' + class_name + '& other) const { return !__eq__(other);}\n')

    override = KEYWORD_OVERRIDE if has_single_superclass(class_def) else ''
    f.write('  ' + RET_TYPE_TO_STR + ' ' + DECL_TO_STR_W_DEFAULT + ' ' + override + ';\n\n')
    
    if decls:
        f.write(decls + '\n\n')
    
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
        # must be set otherwise we get duplicates
        inner_types.add( get_first_inner_type(t) )

    types_no_ptrs = [ t for t in inner_types if not is_yaml_ptr_type(t) ] 
    sorted_types_no_enums = [ t for t in types_no_ptrs if t not in g_enums ]
    sorted_types_no_enums.sort()
    
    for t in sorted_types_no_enums:
        f.write('#include "' + get_api_class_file_name_w_dir(t, EXT_H) + '"\n')

    if KEY_SUPERCLASSES in class_def:
        scs = class_def[KEY_SUPERCLASSES]
        for sc in scs:
            f.write('#include "' + get_api_class_file_name_w_dir(sc, EXT_H) + '"\n')

def write_forward_decls(f, class_def, class_name):
    # first we need to collect all types that we will need
    types = get_all_used_compound_types_no_decorations(class_def)
    
    if class_name not in types and class_def[KEY_TYPE] != VALUE_SUBMODULE:
        f.write('class ' + class_name + ';\n')
        
    for t in types:
        f.write('class ' + t + ';\n')
        
    f.write('class PythonExportContext;\n')
        
    f.write('\n')

    
def generate_class_header(class_name, class_def, decls):
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
            elif is_container_class_no_model(class_name):
                f.write(INCLUDE_API_BASE_EXPORT_CLASS_H + '\n')
        
            if has_superclass_other_than_base(class_def):
                f.write('#include "' + get_api_class_file_name_w_dir(class_def[KEY_SUPERCLASS], EXT_H) + '"\n\n')
        
        write_required_includes(f, class_def)
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        write_forward_decls(f, class_def, class_name)
        
        if has_single_superclass(class_def) or is_container_class_no_model(class_name): 
            write_ctor_define(f, class_def, class_name)
        
        if def_type == VALUE_CLASS:
            write_gen_class(f, class_def, class_name, decls)
        
            f.write('class ' + class_name + ';\n')
        elif def_type == VALUE_SUBMODULE:
            assert decls == ''
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
        f.write('  ' + get_underscored(class_name).upper() + CTOR_SUFFIX + '()\n\n')
        if has_single_superclass(class_def):
            f.write('  // defined in generated file\n') 
        f.write('  ' + RET_TYPE_EXPORT_TO_PYTHON + ' ' + DECL_EXPORT_TO_PYTHON + ' ' + KEYWORD_OVERRIDE + ';\n\n')
        f.write('};\n\n')
        
        f.write(NAMESPACES_END + '\n\n')
        f.write('#endif // ' + guard + '\n')
        
        
def write_is_set_check(f, name, is_list):
    f.write('  if (!is_set(' + name + ')) {\n')
    f.write('    throw ValueError("Parameter \'' + name + '\' must be set')
    if is_list:
        f.write(' and the value must not be an empty list')            
    f.write('.");\n')
    f.write('  }\n')
  

def write_check_semantics_implemetation(f, class_name, items):
    f.write(RET_TYPE_CHECK_SEMANTICS + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_CHECK_SEMANTICS + ' {\n') 
    for attr in items:
        if KEY_DEFAULT not in attr:
            write_is_set_check(f, attr[KEY_NAME], is_yaml_list_type(attr[KEY_TYPE]))
    f.write('}\n\n')    
    
        
def write_to_str_implementation(f, class_name, items, based_on_base_superclass):
    f.write(RET_TYPE_TO_STR + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_TO_STR + ' {\n')

    f.write('  std::stringstream ss;\n')
    if based_on_base_superclass:
        f.write('  ss << get_object_name()')
    else:
        f.write('  ss << "' + class_name + '"')

    if items:
        f.write(' << ": " <<\n')

    last_print_nl = False  
    
    num_attrs = len(items) 
    for i in range(num_attrs):
        name = items[i][KEY_NAME]
        t = items[i][KEY_TYPE]
        
        print_nl = False
        starting_nl = '"' if last_print_nl else '"\\n" << ind + "  " << "'
       
        if is_yaml_list_type(t):
            underlying_type = get_first_inner_type(t)
            if is_yaml_ptr_type(underlying_type):
                f.write('      ' + starting_nl + name + '=" << ' + VEC_PTR_TO_STR + '(' + name + ', all_details, ind + "  ")')
                print_nl = True
            else:
                f.write('      "' + name + '=" << ')
                f.write(VEC_NONPTR_TO_STR + '(' + name + ', all_details, ind + "  ")')
        elif not is_base_yaml_type(t) and t not in g_enums:
            if is_yaml_ptr_type(t):
                f.write('      ' + starting_nl + name + '=" << "(" << ((' + name + ' != nullptr) ? ' + name + '->to_str(all_details, ind + "  ") : "null" ) << ")"')
            else:
                f.write('      ' + starting_nl + name + '=" << "(" << ' + name + '.to_str(all_details, ind + "  ") << ")"')
                
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
    
    f.write('}\n\n')                


def write_vec_export(f, return_only_name, gen_class_name, item):
    
    item_name = item[KEY_NAME];
    PARENT_NAME = 'parent_name'
    
    name_w_args = \
        EXPORT_VEC_PREFIX + item_name + '(' + EXPORT_TO_PYTHON_ARGS + ', const std::string& ' + PARENT_NAME + ')'
    decl = '  virtual ' + RET_TYPE_EXPORT_TO_PYTHON + ' ' + name_w_args + ';\n'
    f.write(RET_TYPE_EXPORT_TO_PYTHON + " " + gen_class_name + "::" + name_w_args + " {\n")
    
    if return_only_name:
        f.write('  // prints vector into \'out\' and returns name of the generated list\n')
    else:
        f.write('  // does not print the array itself to \'out\' and returns the whole list\n')
    
    f.write('  std::stringstream ss;\n')
    out = '  ss << '
    if return_only_name:
        f.write('  std::string ' + EXPORTED_NAME + ';\n' +
                '  if (' + PARENT_NAME + ' != ""){\n' +
                '    ' + EXPORTED_NAME + ' = ' + PARENT_NAME + '+ "_' + item_name + '";\n' +
                '  }\n' +
                '  else {\n' +
                '    ' + EXPORTED_NAME + ' = "' + item_name + '";\n' +
                '  }\n\n')
        f.write(out + EXPORTED_NAME + ' << " = [\\n";\n')
        spacing = '    '
    else:
        f.write(out + '"[";\n')
        spacing = ' ' 

    type = item[KEY_TYPE]
    assert is_yaml_list_type(type)
    
    inner_type = get_inner_list_type(type)

    f.write(
        '  for (size_t i = 0; i < ' + item_name + '.size(); i++) {\n' +
        '    const auto& item = ' + item_name + '[i];\n' +
        '    if (i == 0) {\n' +
        '    ' + out + '"' + spacing + '";\n' +
        '    }\n' +
        '    else if (i % 16 == 0) {\n' +
        '    ' + out + '"\\n  ";\n' +
        '    }\n'
    )
    
    if is_yaml_list_type(inner_type):
        
        inner2_type = get_inner_list_type(inner_type)
        

        # special case for 2D arrays, they hold int or float
        f.write(
            '  ' + out + '"[";\n' +
            '    for (const auto& value: item) {\n'
        )

        if inner2_type == YAML_TYPE_FLOAT:
            f.write('    ' + out + F_TO_STR + '(value) << ", ";\n')
        else:
            f.write('    ' + out + 'value << ", ";\n') 
            
        f.write(
            '    }\n' +
            '  ' + out + '"], ";\n'            
        )
    elif not is_base_yaml_type(inner_type) and inner_type not in g_enums:
        # array of API objects
        f.write(
            '    if (!item->skip_python_export()) {\n' +
            '      std::string name = item->export_to_python(out, ctx);\n' + 
            '    ' + out + 'name << ", ";\n' +
            '    }\n'
        )
    elif inner_type == YAML_TYPE_FLOAT:
        f.write('  ' + out + F_TO_STR + '(item) << ", ";\n')
    elif inner_type == YAML_TYPE_STR:
        f.write('  ' + out + '"\'" << item << "\', ";\n')
    else:
        # array of simple type
        f.write('  ' + out + 'item << ", ";\n')

    f.write('  }\n')
        
    if return_only_name:
        f.write(out + '"\\n]\\n\\n";\n')
        f.write('  out << ss.str();\n')
        f.write('  return ' + EXPORTED_NAME + ';\n')
    else:
        f.write(out + '"]";\n')
        f.write('  return ss.str();\n')
        

    f.write('}\n\n')
    return decl 


def write_export_to_python_implementation(f, class_name, class_def):
    items = class_def[KEY_ITEMS]
    
    # declaration
    if KEY_SUPERCLASS in class_def and class_def[KEY_SUPERCLASS] == BASE_DATA_CLASS:
        override = ' ' + KEYWORD_OVERRIDE
        virtual = ''
    else:
        override = ''
        virtual = KEYWORD_VIRTUAL + ' ' 

    gen_class_name = GEN_CLASS_PREFIX + class_name
    decl_main = '  ' + virtual + RET_TYPE_EXPORT_TO_PYTHON + ' ' + DECL_EXPORT_TO_PYTHON + override + ';\n'
    f.write(RET_TYPE_EXPORT_TO_PYTHON + " " + gen_class_name + "::" + DECL_EXPORT_TO_PYTHON + " {\n")
    
    export_vecs_to_define = []
    
    if class_name == CLASS_NAME_MODEL:
        f.write('# if 0 // not to be used\n')
    
    name_underscored = get_underscored(class_name)
    if not is_container_class(class_name):
        f.write(
            '  if (!export_even_if_already_exported() && ' + CTX + '.already_exported(this)) {\n' +
            '    return ' + CTX + '.get_exported_name(this);\n' +
            '  }\n'
        )
        
        # name generation 
        f.write('  std::string ' + EXPORTED_NAME + ' = ') 
        if not has_item_w_name(items, ATTR_NAME_NAME):
            f.write('"' + name_underscored + '_" + ' + 
                    'std::to_string(' + CTX + '.postinc_counter("' + name_underscored + '"));\n')
        else:
            f.write('std::string("' + name_underscored + '") + "_" + ' + 
                    '(is_set(name) ? ' + 
                        'fix_id(name) : ' +
                        'std::to_string(' + CTX + '.postinc_counter("' + name_underscored + '")));\n')
            
        f.write(
            '  if (!export_even_if_already_exported()) {\n' +
            '    ' + CTX + '.add_exported(this, ' + EXPORTED_NAME + ');\n'
            '  }\n\n'
        )
    else:
        f.write('  std::string ' + EXPORTED_NAME + ' = "' + name_underscored + '";\n\n') 
        
    STR_EXPORT = 'str_export'
    f.write('  bool ' + STR_EXPORT + ' = export_as_string_without_newlines();\n')
    NL = 'nl';
    f.write('  std::string ' + NL + ' = "";\n')  
    IND = 'ind';
    f.write('  std::string ' + IND + ' = " ";\n')  
    
    f.write('  std::stringstream ss;\n')
    out = '  ss << '
    
    f.write('  if (!' + STR_EXPORT + ') {\n')
    f.write('    ' + NL + ' = "\\n";\n')
    f.write('    ' + IND + ' = "    ";\n')
    f.write('  ' + out + EXPORTED_NAME + ' << " = ";\n')
    f.write('  }\n')
    
    f.write(out + '"' + M_DOT + class_name + '(" << ' + NL + ';\n')
    processed_items = set()

    # we would like to export inherited fields first
    sorted_items = [] 
    for item in items:
        if KEY_INHERITED in item:
            sorted_items.append(item)
    for item in items:
        if KEY_INHERITED not in item:
            sorted_items.append(item)
    
    for item in sorted_items:
        
        type = item[KEY_TYPE]
        name = item[KEY_NAME]
        
        if name in processed_items:
            continue
        processed_items.add(name)
        
        export_condition = ''
        if is_yaml_list_type(type):
            export_condition = ' && !skip_vectors_export()'

        if KEY_DEFAULT in item:
            if is_yaml_ptr_type(type):
                f.write('  if (is_set(' + name + ')) {\n')
            else:
                f.write('  if (' + name + ' != ' + get_default_or_unset_value(item) + export_condition + ') {\n')
            f.write('  ')  
            
        f.write(out + IND + ' << "' + name + ' = " << ')
                    
        if is_yaml_list_type(type):
            f.write(EXPORT_VEC_PREFIX + name + '(out, ' + CTX + ', ' + EXPORTED_NAME + ') << "," << ' + NL + ';\n')
            export_vecs_to_define.append(item)
        elif is_yaml_ptr_type(type):
            f.write(name + '->export_to_python(out, ' + CTX + ') << "," << ' + NL + ';\n')
        elif not is_base_yaml_type(type) and type not in g_enums:
            f.write(name + '.export_to_python(out, ' + CTX + ') << "," << ' + NL + ';\n')
        elif type == YAML_TYPE_STR:
            f.write('"\'" << ' + name + ' << "\'" << "," << ' + NL + ';\n')
        elif type == YAML_TYPE_VEC3:
            # also other ivec*/vec2 types should be handled like this but it was not needed so far
            f.write('"' + M_DOT + CPP_TYPE_VEC3 + '(" << ' + 
                        F_TO_STR + '(' + name + '.x) << ", " << ' +
                        F_TO_STR + '(' + name + '.y) << ", " << ' +
                        F_TO_STR + '(' + name + '.z)' +
                        '<< ")," << ' + NL + ';\n')
        elif type == YAML_TYPE_VEC2:
            # also other ivec*/vec2 types should be handled like this but it was not needed so far
            f.write('"' + M_DOT + CPP_TYPE_VEC2 + '(" << ' + 
                        F_TO_STR + '(' + name + '.u) << ", " << ' +
                        F_TO_STR + '(' + name + '.v)' +
                        '<< ")," << ' + NL + ';\n')
        elif type == YAML_TYPE_FLOAT:
            f.write(F_TO_STR + '(' + name + ') << "," << ' + NL + ';\n')
        else:
            # using operators << for enums
            f.write(name + ' << "," << ' + NL + ';\n')
            
        if KEY_DEFAULT in item:
            f.write('  }\n')
            
    f.write(out + '")" << ' + NL + ' << ' + NL + ';\n')
    f.write('  if (!' + STR_EXPORT + ') {\n')
    f.write('    out << ss.str();\n')
    f.write('    return ' + EXPORTED_NAME + ';\n')
    f.write('  }\n')
    f.write('  else {\n')
    f.write('    return ss.str();\n')
    f.write('  }\n')
    
    if class_name == CLASS_NAME_MODEL:
        f.write('#else // # if 0\n')
        f.write('  assert(false);\n')
        f.write('  return "";\n')
        f.write('#endif\n')
        
    f.write('}\n\n')
    
    # also define export vec methods
    decls_export_vec = ''
    for item in export_vecs_to_define:
        decls_export_vec += write_vec_export(f, is_container_class(class_name), gen_class_name, item)
    
    # return string that contains declaration of export methods, this will be later used in 
    # .h file generation
    return decl_main + decls_export_vec


def write_operator_equal_body(f, class_name, class_def, skip_arrays_and_name=False):
    items = class_def[KEY_ITEMS]

    f.write('  return\n') 

    if not items:
        f.write('true ;\n')
   
    num_attrs = len(items) 
    for i in range(num_attrs):
        name = items[i][KEY_NAME]
        t = items[i][KEY_TYPE]
        if is_yaml_ptr_type(t):
            
            f.write(
                '    (\n' 
                '      (is_set(' + name + ')) ?\n'  
                '        (is_set(other.' + name + ') ?\n'
                '          (' + name + '->__eq__(*other.' + name + ')) : \n'
                '          false\n' 
                '        ) :\n'
                '        (is_set(other.' + name + ') ?\n'
                '          false :\n'
                '          true\n'
                '        )\n'
                '     ) '
            )
            
        elif is_yaml_list_type(t):
            if not skip_arrays_and_name:
                if is_yaml_ptr_type(get_inner_list_type(t)):
                    f.write('    vec_ptr_eq(' + name + ', other.' + name + ')')
                else:
                    f.write('    ' + name + ' == other.' + name)
            else:
                f.write('    true /*' + name + '*/')
        else:
            if skip_arrays_and_name and name == KEY_NAME:
                f.write('    (ignore_name || ' + name + ' == other.' + name + ')')
            else:
                f.write('    ' + name + ' == other.' + name)
            
        if i != num_attrs - 1:
            f.write(' &&')
        else: 
            f.write(';')
        f.write('\n')
       
       
def write_copy_implementation(f, class_name, class_def, deepcopy):

    if deepcopy:
        func_name = get_deepcopy_function_name(class_name)
        func_args = 'py::dict'
    else:
        func_name = get_copy_function_name(class_name)
        func_args = ''
    
    f.write(get_as_shared_ptr(class_name) + ' ' + GEN_CLASS_PREFIX + class_name + '::' + func_name + '(' + func_args + ') const {\n')
    
    # for some reason res(DefaultCtorArgType()) is not accepted by gcc...
    f.write('  ' + SHARED_PTR + '<' + class_name + '> res = ' + MAKE_SHARED + '<' + class_name + '>(' + DEFAULT_CTOR_ARG_TYPE +'());\n')

    items = class_def[KEY_ITEMS]

    if has_single_superclass(class_def):
        f.write('  res->' + CLASS_NAME_ATTR + ' = ' + CLASS_NAME_ATTR + ';\n')
    
    # TODO - deepcopy - for 
    for item in items:
        name = item[KEY_NAME]
        
        if not deepcopy:
            # simply use assign operator
            f.write('  res->' + name + ' = ' + name + ';\n')
        else:
            # pointers must be deepcopied
            # also vectors containing pointers, cannot use aux functions for copying vectors 
            # because name of the deepcpy function is different for each class
            t = item[KEY_TYPE]
            
            # check for 2D lists
            innermost_list_ptr_type = None
            if is_yaml_list_type(t):
                inner = get_inner_list_type(t)
                if is_yaml_ptr_type(inner):
                    innermost_list_ptr_type = inner
                elif is_yaml_list_type(inner):
                    inner2 = get_inner_list_type(inner)
                    if is_yaml_ptr_type(inner2):
                        assert False, "2D lists containing pointers are not supported for deepcopy yet"
            
            if innermost_list_ptr_type:
                base_t = remove_ptr_mark(innermost_list_ptr_type)
                f.write('  for (const auto& item: ' + name + ') {\n')
                f.write('    res->' + name + '.push_back((' + IS_SET + '(item)) ? ' + 
                             'item->' + get_deepcopy_function_name(base_t) + '() : '
                             'nullptr);\n')
                f.write('  }\n')
            elif is_yaml_dict_type(t):
                assert False, "Dict type is not supported for deepcopy yet"
            elif is_yaml_ptr_type(t):
                base_t = remove_ptr_mark(t)
                f.write('  res->' + name + ' = ' + IS_SET + '(' + name + ') ? ' + 
                             name + '->' + get_deepcopy_function_name(base_t) + '() : '
                             'nullptr;\n')
            else:
                f.write('  res->' + name + ' = ' + name + ';\n')

    f.write('\n')
    f.write('  return res;\n')
    f.write('}\n\n')
        
        
def write_operator_equal_implementation(f, class_name, class_def):
    # originally was operator== used, but this causes too cryptic compilation errors, 
    # so it was better to use python-style naming,
    # also the __eq__ is only defined for the 'Gen' classes, so defining operator ==
    # might have been more confusing 
                
    gen_class_name = GEN_CLASS_PREFIX + class_name
    f.write('bool ' + gen_class_name + '::__eq__(const ' + class_name + '& other) const {\n')
    write_operator_equal_body(f, class_name, class_def)
    f.write('}\n\n')    
    
    f.write('bool ' + gen_class_name + '::eq_nonarray_attributes(const ' + class_name + '& other, const bool ignore_name) const {\n')
    write_operator_equal_body(f, class_name, class_def, skip_arrays_and_name=True)
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


def create_doc_str(yaml_doc, w_quotes=True):
    if not yaml_doc:
        return ""
    
    res = yaml_doc.replace('"', '\\"')
    res = res.replace('\n', '\\n')
    res = res.replace('\\:', ':')
    if w_quotes:
        return '"' + res + '"'
    else:
        return res


def write_pybind11_method_bindings(f, class_name, method, class_def):
    assert KEY_NAME in method
    name = rename_cpp_reserved_id(method[KEY_NAME])
    
    full_method_name = '&' + class_name + '::' + name
    if is_overloaded(method, class_def):
        # overloaded method must be extra decorated with argument types for pybind11
        full_method_name = get_method_overload_cast(method) + '(' + full_method_name + ')' 
    
    f.write('      .def("' + name + '", ' + full_method_name)  
    
    params_doc = ''
    if KEY_PARAMS in method:
        params = method[KEY_PARAMS]
        params_cnt = len(params)
        for i in range(params_cnt):
            f.write(', ')
            p = params[i]
            assert KEY_NAME in p
            param_name = p[KEY_NAME]
            f.write('py::arg("' + param_name + '")')
            if KEY_DEFAULT in p:
                f.write(' = ' + get_default_or_unset_value(p))
            
            params_doc += '- ' + param_name
            if KEY_DOC in p:
                params_doc += ': ' + create_doc_str(p[KEY_DOC], False) + '\\n'
            params_doc += '\\n'
            
            
    if KEY_DOC in method or params_doc:
        
        f.write(', "') 
        if KEY_DOC in method:
            f.write(create_doc_str(method[KEY_DOC], False))
            if params_doc:
                params_doc = '\\n' + params_doc
        f.write(params_doc + '"')
                
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
                SHARED_PTR + '<' + class_name + '>>(m, "' + class_name + '"')
        if KEY_DOC in class_def:
            f.write(', ' + create_doc_str(class_def[KEY_DOC]))                
                 
        f.write(')\n')
        f.write('      .def(\n')
        f.write('          py::init<\n')
    
        num_items = len(items)
        if ctor_has_args:
            if has_single_superclass(class_def) or is_container_class_no_model(class_name):
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
            if has_single_superclass(class_def) or is_container_class_no_model(class_name):
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
        
        f.write('      .def("__copy__", &' + class_name + '::' + get_copy_function_name(class_name) + ')\n')
        f.write('      .def("__deepcopy__", &' + class_name + '::' + get_deepcopy_function_name(class_name) + ', py::arg("memo"))\n')
        f.write('      .def("__str__", &' + class_name + '::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))\n')
        # keeping the default __repr__ implementation for better error messages 
        #f.write('      .def("__repr__", &' + class_name + '::to_str, py::arg("ind") = std::string(""))\n')
        f.write('      .def("__eq__", &' + class_name + '::__eq__, py::arg("other"))\n')
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
                f.write('      .def_property("' + name + '", &' + class_name + '::get_' + name + ', &' + class_name + '::set_' + name)
                
                if is_yaml_list_type(items[i][KEY_TYPE]):
                    f.write(', py::return_value_policy::reference') 
                
                if KEY_DOC in items[i]:
                    f.write(', ' + create_doc_str(items[i][KEY_DOC]))
        
                f.write(')\n')
                
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
    # returns a list of methods to be declared 
    decls = ''
    with open(get_gen_class_file_name_w_dir(class_name, EXT_CPP), 'w') as f:
        def_type = class_def[KEY_TYPE]
        
        f.write(COPYRIGHT)
        write_generated_notice(f)
        
        f.write('#include <sstream>\n')
        f.write('#include "api/pybind11_stl_include.h"\n')
        f.write(INCLUDE_API_PYTHON_EXPORT_UTILS_H + '\n')

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
                write_set_initialized_implemetation(f, class_name, items)
                write_set_all_default_or_unset(f, class_name, class_def)
    
            write_copy_implementation(f, class_name, class_def, False)
            write_copy_implementation(f, class_name, class_def, True)
            write_operator_equal_implementation(f, class_name, class_def)
            write_to_str_implementation(f, class_name, items, has_single_superclass(class_def))
        
        write_pybind11_bindings(f, class_name, class_def)

        # we do not need export for introspected classes
        if is_container_class(class_name) or \
            (KEY_SUPERCLASS in class_def and class_def[KEY_SUPERCLASS] != BASE_INTROSPECTION_CLASS):
            decls += write_export_to_python_implementation(f, class_name, class_def);
        
        f.write(NAMESPACES_END + '\n\n')
        
    return decls 

# class Model provides the same methods as Subsystem and Instantiation, 
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
    
    decls = generate_class_implementation_and_bindings(class_name, class_def_w_inheritances)
    generate_class_header(class_name, class_def_w_inheritances, decls)
    
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

    
    f.write('\nstatic inline std::ostream& operator << (std::ostream& out, const ' + name + ' v) {\n')
    f.write('  switch (v) {\n');
    for i in range(num):
        v = values[i]
        # dumping in Python style '.' that will be most probably more common as API  
        #f.write('    case ' + name + '::' + v[KEY_NAME] + ': out << "' + name + '.' + v[KEY_NAME] + ' (' + str(v[KEY_VALUE]) + ')"; break;\n')
        
        # dumping as Python code  
        f.write('    case ' + name + '::' + v[KEY_NAME] + ': out << "' + M_DOT + name + '.' + v[KEY_NAME] + '"; break;\n')
        
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
        f.write('#include <climits>\n')
        f.write('#include <cfloat>\n')

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
    num = len(values)
    
    f.write('  py::enum_<' + name + '>(m, "' + name + '", py::arithmetic()')

    # Pybind11 does nto allow help for enum members so we put all that information into
    # the header 
    members_doc = ''
    for i in range(num):
        v = values[i]
        assert KEY_NAME in v
        members_doc += '- ' + v[KEY_NAME]
        if KEY_DOC in v:
            members_doc += ': ' + create_doc_str(v[KEY_DOC], False)
        else:
            members_doc += '\\n'
        members_doc += '\\n'

    if KEY_DOC in enum_def or members_doc:
        f.write(', "')
        if KEY_DOC in enum_def: 
            f.write(create_doc_str(enum_def[KEY_DOC], False) + '\\n')
        f.write(members_doc + '"')

    f.write(')\n')
    
    
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


def capitalize_first(s):
    if s[0].islower():
        return s[0].upper() + s[1:]
    else:
        return s


def generate_vector_bindings(data_classes):
    # collect all attributes that use List
    list_types = set()
    used_types = set()
    for key, value in data_classes.items():
        
        if key != KEY_CONSTANTS and key != KEY_ENUMS:
            # items 
            if KEY_ITEMS in value:
                for item in value[KEY_ITEMS]:
                    if is_yaml_list_type(item[KEY_TYPE]):
                        t = item[KEY_TYPE]
                        list_types.add(t)
                        
                        inner = get_inner_list_type(t)
                        if is_yaml_list_type(inner):
                            inner = get_inner_list_type(inner)
                        used_types.add(inner)
                        
    sorted_list_types = sorted(list_types)
    sorted_used_types = sorted(used_types)
    
    # generate gen_make_opaque.h
    with open(os.path.join(TARGET_DIRECTORY, GEN_VECTORS_MAKE_OPAQUE_H), 'w') as f:
        f.write(COPYRIGHT)
        guard = 'GEN_VECTORS_MAKE_OPAQUE_H'
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        f.write('#include <vector>\n')
        f.write('#include <memory>\n')
        f.write('#include "pybind11/include/pybind11/pybind11.h"\n') 
        f.write('#include "defines.h"\n\n')
        f.write(NAMESPACES_BEGIN + '\n\n')
        
        for t in sorted_used_types:
            if is_base_yaml_type(t):
                continue
            assert is_yaml_ptr_type(t)
            f.write('class ' + t[:-1] + ';\n')
        
        f.write('\n' + NAMESPACES_END + '\n\n')
        
        for t in sorted_list_types:
            f.write(PYBIND11_MAKE_OPAQUE + '(' + yaml_type_to_cpp_type(t, True) + ');\n')
        
        f.write('\n#endif // ' + guard + '\n')
             
    # generate gen_bind_vector.cpp
    with open(os.path.join(TARGET_DIRECTORY, GEN_VECTORS_BIND_CPP), 'w') as f:
        f.write(COPYRIGHT + '\n')
        f.write('#include "api/pybind11_stl_include.h"\n')
        f.write('#include "pybind11/include/pybind11/stl_bind.h"\n')
        f.write('#include "pybind11/include/pybind11/numpy.h"\n\n')
        
        f.write('namespace py = pybind11;\n\n')
        
        for t in sorted_used_types:
            if is_base_yaml_type(t):
                continue
            assert is_yaml_ptr_type(t)
            f.write('#include "' +  get_api_class_file_name_w_dir(t[:-1], EXT_H) + '"\n')
        f.write('\n')
        
        f.write(NAMESPACES_BEGIN + '\n\n')

        f.write('void ' + GEN_VECTORS_BIND + '{\n')
        
        ind = '  '
        for t in sorted_list_types:
            cpp_t = yaml_type_to_cpp_type(t, True)
            
            name = 'Vector'
            inner = get_inner_list_type(t)
            if is_yaml_list_type(inner):
                inner2 = get_inner_list_type(inner)
                name += name + capitalize_first(inner2)
            else:
                name += capitalize_first(inner)
                
            if name[-1] == '*':
                name = name[:-1]
                            
            f.write(ind + PY_BIND_VECTOR + '<' + cpp_t + '>(m,"' + name + '");\n')
            f.write(ind + PY_IMPLICITLY_CONVERTIBLE + '<py::list, ' + cpp_t + '>();\n')
            f.write(ind + PY_IMPLICITLY_CONVERTIBLE + '<py::tuple, ' + cpp_t + '>();\n')
        
            # special case for numpy conversions
            if name == 'VectorFloat':
                f.write(ind + PY_IMPLICITLY_CONVERTIBLE + '<py::array_t<double>, ' + cpp_t + '>();\n')
            elif name == 'VectorInt':
                f.write(ind + PY_IMPLICITLY_CONVERTIBLE + '<py::array_t<int>, ' + cpp_t + '>();\n')
            f.write('\n')
                
        
        f.write('}\n')
        f.write('\n' + NAMESPACES_END + '\n\n')
        

def generate_data_classes(data_classes):
    generate_constants_and_enums(
        data_classes[KEY_CONSTANTS] if KEY_CONSTANTS in data_classes else [],
        data_classes[KEY_ENUMS] if KEY_ENUMS in data_classes else [])

    set_global_enums(data_classes)
    
    # std::vector/list needs special handling because we need it to be opaque 
    # so that the user can modify the internal value
    generate_vector_bindings(data_classes)

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
        
        f.write('INT32_MAX = 2147483647 # do not use this constant in your code\n\n')
        f.write('FLT_MAX = 3.40282346638528859812e+38 # do not use this constant in your code\n\n')
        
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
            'class Vec2():\n'
            '    def __init__(self, u : float = 0, v : float = 0):\n'
            '        self.u = u\n'
            '        self.v = v\n'
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
        
        
def indent_and_fix_rst_chars(text, indent):
    return text.replace('\n', '\n' + indent).replace('*', '\\*').replace('@', '\\@')
            

def load_and_generate_data_classes():
    # create work directory
    work_dir = os.path.join('..', 'work')
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    
    data_classes = {}
    
    for cat in CATEGORIES:
        input_file = cat + '.yaml'
        with open(input_file) as file:
            # The FullLoader parameter handles the conversion from YAML
            # scalar values to Python the dictionary format
            loaded_data = yaml.load(file, Loader=yaml.FullLoader)
            if loaded_data:
                # set information about the original file name
                for key,value in loaded_data.items():
                    if key != KEY_CONSTANTS and key != KEY_ENUMS:
                        value[KEY_CATEGORY] = cat
                 
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
    doc.generate_documentation(data_classes)

if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == '-v':
        VERBOSE = True
        
    load_and_generate_data_classes()
    
    