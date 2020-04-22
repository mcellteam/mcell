#!/usr/bin/env python3
"""
Copyright (C) 2020 by
The Salk Institute for Biological Studies and
Pittsburgh Supercomputing Center, Carnegie Mellon University

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

import sys
import os
import yaml
import re
from datetime import datetime
from copy import copy

DATA_CLASSES_FILE = 'data_classes.yaml'


COPYRIGHT = \
"""/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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
API_DIRECTORY = os.path.join('..', 'api')

KEY_ITEMS = 'items'
KEY_NAME = 'name'
ATTR_NAME_NAME = 'name' # attrribute with name 'name' is already defined in BaseDataClass
KEY_TYPE = 'type'
KEY_VALUE = 'value'
KEY_DEFAULT = 'default'

KEY_SUPERCLASS = 'superclass'
KEY_SUPERCLASSES = 'superclasses'

KEY_CONSTANTS = 'constants'

KEY_METHODS = 'methods'
KEY_PARAMS = 'params'
KEY_RETURN_TYPE = 'return_type'

KEY_INHERITED = 'inherited' # used only internally, not in input YAML

YAML_TYPE_FLOAT = 'float'
YAML_TYPE_STR = 'str'
YAML_TYPE_INT = 'int'
YAML_TYPE_LONG = 'long'
YAML_TYPE_VEC2 = 'Vec2'
YAML_TYPE_VEC3 = 'Vec3'
YAML_TYPE_LIST = 'List'

CPP_TYPE_FLOAT = 'float_t'
CPP_TYPE_STR = 'std::string'
CPP_TYPE_INT = 'int'
CPP_TYPE_LONG = 'long'
CPP_TYPE_VEC2 = 'Vec2'
CPP_TYPE_VEC3 = 'Vec3'
CPP_VECTOR_TYPE = 'std::vector'

CPP_REFERENCE_TYPES = [CPP_TYPE_STR, CPP_TYPE_VEC2, CPP_TYPE_VEC3]

UNSET_VALUE = 'unset'
EMPTY_ARRAY = 'empty'

UNSET_VALUE_FLOAT = 'FLT_UNSET'
UNSET_VALUE_STR = 'STR_UNSET'
UNSET_VALUE_INT = 'INT_UNSET'
UNSET_VALUE_INT = 'LONG_UNSET'
UNSET_VALUE_VEC2 = 'VEC2_UNSET'
UNSET_VALUE_VEC3 = 'VEC3_UNSET'
UNSET_VALUE_PTR = 'nullptr'


GEN_PREFIX = 'gen_'
GUARD_PREFIX = 'API_GEN_'
GUARD_SUFFIX = '_H'
CTOR_SUFFIX = '_CTOR'
EXT_CPP = 'cpp'
EXT_H = 'h'

GEN_CLASS_PREFIX = 'Gen'
BASE_DATA_CLASS = 'BaseDataClass'

RET_TYPE_CHECK_SEMANTICS = 'SemRes'
DECL_CHECK_SEMANTICS = 'check_semantics(std::ostream& out) const'
DECL_DEFINE_PYBINDIND_CONSTANTS = 'void define_pybinding_constants(py::module& m)'
RET_TYPE_TO_STR = 'std::string'
SHARED_PTR = 'std::shared_ptr'
DECL_TO_STR = 'to_str() const'
KEYWORD_OVERRIDE = 'override'
  
CLASS_NAME_ATTR = 'class_name'

GEN_CONSTANTS_H = 'gen_constants.h'
GEN_CONSTANTS_CPP = 'gen_constants.cpp'

INCLUDE_API_MCELL_H = '#include "../api/mcell.h"'
INCLUDE_API_COMMON_H = '#include "../api/common.h"'
NAMESPACES_BEGIN = 'namespace MCell {\nnamespace API {'
NAMESPACES_END = '} // namespace API\n} // namespace MCell'

VEC_NONPTR_TO_STR = 'vec_nonptr_to_str'
VEC_PTR_TO_STR = 'vec_ptr_to_str'

def get_underscored(class_name):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', class_name).lower()

def get_header_guard_name(class_name):
    return GUARD_PREFIX + get_underscored(class_name).upper() + GUARD_SUFFIX

def get_gen_class_file_name(class_name, extension):
    return os.path.join(GEN_PREFIX + get_underscored(class_name) + '.' + extension)

def get_gen_class_file_name_w_dir(class_name, extension):
    return os.path.join(TARGET_DIRECTORY, get_gen_class_file_name(class_name, extension))

def get_api_class_file_name(class_name, extension):
    return os.path.join(API_DIRECTORY, get_underscored(class_name) + '.' + extension)

def is_list(t):
    return t.startswith(YAML_TYPE_LIST)


# rename inner to underlying?
def get_inner_list_type(t):
    if is_list(t):
        return t[len(YAML_TYPE_LIST)+1:-1]
    else:
        return t


def is_base_yaml_type(t):
    return \
        t == YAML_TYPE_FLOAT or t == YAML_TYPE_STR or t == YAML_TYPE_INT or t == YAML_TYPE_LONG or \
        t == YAML_TYPE_VEC2 or t == YAML_TYPE_VEC3 or \
        (is_list(t) and is_base_yaml_type(get_inner_list_type(t)))


def is_yaml_ptr_type(t):
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
    elif t == YAML_TYPE_VEC2:
        return CPP_TYPE_VEC2
    elif t == YAML_TYPE_VEC3:
        return CPP_TYPE_VEC3
    elif is_list(t):
        assert len(t) > 7
        inner_type = yaml_type_to_cpp_type(get_inner_list_type(t))
        return CPP_VECTOR_TYPE + '<' + inner_type + '>'
    else:
        if is_yaml_ptr_type(t):
            return SHARED_PTR + '<' + t[0:-1] + '>' 
        else:
            return t # standard ttype
    
    
def yaml_type_to_pybind_type(t):
    assert len(t) >= 1
    if t == YAML_TYPE_FLOAT:
        return CPP_TYPE_FLOAT
    elif t == YAML_TYPE_STR:
        return YAML_TYPE_STR
    elif t == YAML_TYPE_INT:
        return CPP_TYPE_INT + '_'
    elif t == YAML_TYPE_LONG:
        return CPP_TYPE_INT + '_'
    else:
        assert False, "Unsupported constant type " + t 
    
    
def get_type_as_ref_param(attr):
    assert KEY_TYPE in attr
    
    yaml_type = attr[KEY_TYPE]
    cpp_type = yaml_type_to_cpp_type(yaml_type)
    
    res = cpp_type
    if cpp_type in CPP_REFERENCE_TYPES:
        res += '&'
    
    return res


def get_default_or_unset_value(attr):
    if KEY_DEFAULT in attr:
        default_value = attr[KEY_DEFAULT]
        if default_value != UNSET_VALUE and default_value != EMPTY_ARRAY:
            return str(default_value) 
    
    assert KEY_TYPE in attr
    t = attr[KEY_TYPE]
    if t == YAML_TYPE_FLOAT:
        return UNSET_VALUE_FLOAT
    elif t == YAML_TYPE_STR:
        return UNSET_VALUE_STR
    elif t == YAML_TYPE_INT:
        return UNSET_VALUE_INT
    elif t == YAML_TYPE_LONG:
        return UNSET_VALUE_LONG
    elif t == YAML_TYPE_VEC2:
        return UNSET_VALUE_VEC2
    elif t == YAML_TYPE_VEC3:
        return UNSET_VALUE_VEC3
    elif is_list(t):
        return yaml_type_to_cpp_type(t) + '()'
    else:
        return UNSET_VALUE_PTR


def is_cpp_ptr_type(t):
    return t.startswith(SHARED_PTR)

    
def is_cpp_ref_type(t):
    return t in CPP_REFERENCE_TYPES


def write_generated_notice(f, input_file_name):
    now = datetime.now()
    # date printing is disabled during the development phase to minimize changes in files 
    #date_time = now.strftime("%m/%d/%Y, %H:%M")
    #f.write('// This file was generated automatically on ' + date_time + ' from ' + '\'' + input_file_name + '\'\n\n')

def write_ctor_body(f, class_def, class_name, append_backslash, indent, only_inherited):
    items = class_def[KEY_ITEMS] if KEY_ITEMS in class_def else []
    backshlash = '\\' if append_backslash else ''
    
    f.write(indent + class_name + '( ' + backshlash + '\n')
     
    inherited_items = [ attr for attr in items if is_inherited(attr) ]

    if only_inherited:
        # override also the original items list
        items = inherited_items 
     
    # ctor parameters    
    num_items = len(items)
    for i in range(num_items):
        attr = items[i]
        
        assert KEY_NAME in attr
        name = attr[KEY_NAME]

        const_spec = 'const ' if not is_yaml_ptr_type(attr[KEY_TYPE]) else ''
        f.write(indent + '    ' + const_spec + get_type_as_ref_param(attr) + ' ' + name + '_')
        
        if KEY_DEFAULT in attr:
            f.write(' = ' + get_default_or_unset_value(attr))
        
        if i != num_items - 1:
            f.write(',')
        f.write(' ' + backshlash + '\n')

    f.write(indent + ') ')
    
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
    f.write('#define ' + get_underscored(class_name).upper() + CTOR_SUFFIX + '() \\\n')

    write_ctor_body(f, class_def, class_name, append_backslash=True, indent='    ', only_inherited=False)

    f.write('{ \\\n')

    # initialization code
    f.write('      ' + CLASS_NAME_ATTR + ' = "' + class_name + '"; \\\n')
    items = class_def[KEY_ITEMS] if KEY_ITEMS in class_def else []    
    num_items = len(items)
    for i in range(num_items):
        assert KEY_NAME in items[i] 
        attr_name = items[i][KEY_NAME]
        f.write('      ' + attr_name + ' = ' + attr_name + '_; \\\n')
    f.write('    }\n\n')    
    
    
def write_ctor_for_superclass(f, class_def, class_name):
    write_ctor_body(f, class_def, GEN_CLASS_PREFIX + class_name, append_backslash=False, indent='  ', only_inherited=True)
    f.write(' {\n')
    f.write('  }\n')
    

def write_attr_with_get_set(f, attr):
    assert KEY_NAME in attr
    assert KEY_TYPE in attr
    
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
    f.write(method[KEY_NAME] + '(')
    
    if KEY_PARAMS in method:
        params = method[KEY_PARAMS]
        params_cnt = len(params)
        for i in range(params_cnt):
            p = params[i]
            assert KEY_NAME in p
            assert KEY_TYPE in p
            t = get_type_as_ref_param(p)
            const_spec = 'const ' if not is_yaml_ptr_type(p[KEY_TYPE]) else ''
            f.write(const_spec + t + ' ' + p[KEY_NAME])
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
    

def has_superclass(class_def):
    if KEY_SUPERCLASS in class_def:
        return True
    else:
        return False


def has_superclass_other_than_base(class_def):
    # TODO: also deal with superclasses?
    if has_superclass(class_def):
        return class_def[KEY_SUPERCLASS] != BASE_DATA_CLASS
    else:
        return False
    
    
def write_gen_class(f, class_def, class_name):
    f.write('class ' + GEN_CLASS_PREFIX + class_name)
    
    if has_superclass(class_def):
        f.write(': public ' + class_def[KEY_SUPERCLASS])
        
    f.write( ' {\n')
    f.write('public:\n')
    
    if has_superclass_other_than_base(class_def):
        write_ctor_for_superclass(f, class_def, class_name)
    
    if not has_superclass(class_def):
        # generate virtual destructor
        f.write('  virtual ~' + GEN_CLASS_PREFIX + class_name + '() {}\n')
        
    if has_superclass(class_def):
        f.write('  ' + RET_TYPE_CHECK_SEMANTICS + ' ' + DECL_CHECK_SEMANTICS + ' ' + KEYWORD_OVERRIDE + ';\n')
        f.write('  ' + RET_TYPE_TO_STR + ' ' + DECL_TO_STR + ' ' + KEYWORD_OVERRIDE + ';\n\n')
        
    f.write('  // --- attributes ---\n')
    items = class_def[KEY_ITEMS]
    for attr in items:
        if not is_inherited(attr):
            written = write_attr_with_get_set(f, attr)
            if written :
                f.write('\n')

    f.write('  // --- methods ---\n')
    methods = class_def[KEY_METHODS]
    for m in methods:
        if not is_inherited(m):
            write_method_declaration(f, m)
        
    f.write('}; // ' + GEN_CLASS_PREFIX + class_name + '\n\n')


def write_define_binding_decl(f, class_name):
    f.write('py::class_<' + class_name + '> define_pybinding_' + class_name + '(py::module& m)')           


def remove_ptr_mark(t):
    assert len(t) > 1
    if t[-1] == '*':
        return t[0:-1]
    else:
        return t
    

def get_all_used_compound_types(class_def):
    types = set()
    for items in class_def[KEY_ITEMS]:
        assert KEY_TYPE in items
        t = items[KEY_TYPE]
        if not is_base_yaml_type(t):
            types.add( remove_ptr_mark(t) )

    for method in class_def[KEY_METHODS]:
        if KEY_RETURN_TYPE in method:
            t = method[KEY_RETURN_TYPE]
            if not is_base_yaml_type(t):
                types.add( remove_ptr_mark(t) )
            
        if KEY_PARAMS in method:
            for param in method[KEY_PARAMS]:
                assert KEY_TYPE in param
                t = param[KEY_TYPE]
                if not is_base_yaml_type(t):
                    types.add( remove_ptr_mark(t) )
                    
    sorted_types = [ remove_ptr_mark( get_inner_list_type(t) ) for t in list(types) ]
    sorted_types.sort()
    return sorted_types


def write_forward_decls(f, class_def):
    # first we need to collect all types that we will need
    types = get_all_used_compound_types(class_def)
    
    for t in types:
        f.write('class ' + t + ';\n')
        
    if types:
        f.write('\n')
    

def generate_class_header(class_name, class_def, input_file_name):
    with open(get_gen_class_file_name_w_dir(class_name, EXT_H), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f, input_file_name)
        
        guard = get_header_guard_name(class_name);
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        f.write(INCLUDE_API_COMMON_H + '\n')
        
        if has_superclass_other_than_base(class_def):
            f.write('#include "' + get_api_class_file_name(class_def[KEY_SUPERCLASS], EXT_H) + '"\n\n')
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        write_forward_decls(f, class_def)
        
        if has_superclass(class_def): # not sure about this condition
            write_ctor_define(f, class_def, class_name)
        
        write_gen_class(f, class_def, class_name)
        
        f.write('class ' + class_name + ';\n')
        write_define_binding_decl(f, class_name)
        f.write(';\n')
        
        f.write(NAMESPACES_END + '\n\n')
        f.write('#endif // ' + guard + '\n')
    

def write_is_set_check(f, name):
    f.write('  if (!is_set(' + name + ')) {\n')
    f.write('    out << get_object_name() << ": Parameter \'' + name + '\' must be set.\\n";\n')
    f.write('    return SemRes::ERROR;\n')
    f.write('  }\n')
  

def write_check_semantics_implemetation(f, class_name, items):
    f.write(RET_TYPE_CHECK_SEMANTICS + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_CHECK_SEMANTICS + ' {\n') 
    for attr in items:
        if KEY_DEFAULT not in attr:
            write_is_set_check(f, attr[KEY_NAME])
    
    f.write('  return SemRes::OK;\n')
    f.write('}\n\n')    
    
        
def write_to_str_implemetation(f, class_name, items):
    f.write(RET_TYPE_TO_STR + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_TO_STR + ' {\n')
    f.write('  std::stringstream ss;\n')
    f.write('  ss << get_object_name() << ": " <<\n')
    
    num_attrs = len(items) 
    for i in range(num_attrs):
        name = items[i][KEY_NAME]
        f.write('      "' + name + '=" << ')
        
        type = items[i][KEY_TYPE]
        if is_list(type):
            underlying_type = get_inner_list_type(type)
            if is_yaml_ptr_type(underlying_type):
                f.write(VEC_PTR_TO_STR + '(' + name + ')')
            else:
                f.write(VEC_NONPTR_TO_STR + '(' + name + ')')
        elif not is_base_yaml_type(type):
            f.write('"(" << ((' + name + ' != nullptr) ? ' + name + '->to_str() : "null" ) << ")"')
        else:
            f.write(name)

        if i != num_attrs - 1:
            f.write(' << ", " <<\n')
            
    f.write(';\n')
    
    f.write('  return ss.str();\n')
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
    name = method[KEY_NAME]
    
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
                q = '"' if p[KEY_TYPE] == YAML_TYPE_STR else ''
                f.write(' = ' + q + get_default_or_unset_value(p) + q)
                
    f.write(')\n')
    
    
def write_pybind11_bindings(f, class_name, class_def):
    items = class_def[KEY_ITEMS]
    
    write_define_binding_decl(f, class_name)
    f.write(' {\n')
    f.write('  return py::class_<' + class_name + ', ' + SHARED_PTR + '<' + class_name + '>>(m, "' + class_name + '")\n')
    f.write('      .def(\n')
    f.write('          py::init<\n')

    # init operands
    num_items = len(items)
    for i in range(num_items):
        attr = items[i]
        const_spec = 'const ' if not is_yaml_ptr_type(attr[KEY_TYPE]) else ''
        f.write('            ' + const_spec + get_type_as_ref_param(attr))
        if i != num_items - 1:
            f.write(',\n')
    if num_items != 0:
        f.write('\n')
    f.write('          >()\n')
        
    if num_items != 0:
        f.write(',')
    
    # init argument names and default values
    for i in range(num_items):
        attr = items[i]
        name = attr[KEY_NAME]
        f.write('          py::arg("' + name + '")')
        if KEY_DEFAULT in attr:
            f.write(' = ' + get_default_or_unset_value(attr))
        if i != num_items - 1:
            f.write(',\n')
    f.write('\n')            
    f.write('        )\n')            
    
    # common methods
    if has_superclass(class_def):
        f.write('      .def("check_semantics", &' + class_name + '::check_semantics_cerr)\n')
        f.write('      .def("__str__", &' + class_name + '::to_str)\n')
        
    # declared methods
    for m in class_def[KEY_METHODS]:
        write_pybind11_method_bindings(f, class_name, m, class_def)

    # dump needs to be always implemented
    f.write('      .def("dump", &' + class_name + '::dump)\n')
    
    # properties
    for i in range(num_items):
        name = items[i][KEY_NAME]
        f.write('      .def_property("' + name + '", &' + class_name + '::get_' + name + ', &' + class_name + '::set_' + name + ')\n')
    f.write('    ;\n')
    f.write('}\n\n')
    
                
def write_used_classes_includes(f, class_def):
    types = get_all_used_compound_types(class_def)
    for t in types:
        f.write('#include "' + get_api_class_file_name(t, EXT_H) + '"\n')

            
def generate_class_implementation_and_bindings(class_name, class_def, input_file_name):
    with open(get_gen_class_file_name_w_dir(class_name, EXT_CPP), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f, input_file_name)
        
        f.write('#include <sstream>\n')
        f.write('#include <pybind11/stl.h>\n')

        # includes for our class
        f.write('#include "' + get_gen_class_file_name(class_name, EXT_H) + '"\n')  
        f.write('#include "' + get_api_class_file_name(class_name, EXT_H) + '"\n')
        
        # we also need includes for every type that we used
        write_used_classes_includes(f, class_def)
        
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        if has_superclass(class_def):
            items = class_def[KEY_ITEMS]
            write_check_semantics_implemetation(f, class_name, items)
            write_to_str_implemetation(f, class_name, items)
        
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

        if KEY_ITEMS in superclass_def:
            for item in superclass_def[KEY_ITEMS]:
                item_copy = copy(item)
                item_copy[KEY_INHERITED] = True
                res[KEY_ITEMS].append(item_copy) 
        
    return res


def generate_class_files(data_classes, class_name, class_def, input_file_name):
    # we need items and methods to be present 
    if KEY_ITEMS not in class_def:
        class_def[KEY_ITEMS] = []
        
    if KEY_METHODS not in class_def:
        class_def[KEY_METHODS] = []

    class_def_w_inheritances = inherit_from_superclasses(data_classes, class_name, class_def)
    
    generate_class_header(class_name, class_def_w_inheritances, input_file_name)
    generate_class_implementation_and_bindings(class_name, class_def_w_inheritances, input_file_name)
    
    
def write_constant_def(f, constant_def):
    assert KEY_NAME in constant_def
    assert KEY_TYPE in constant_def
    assert KEY_VALUE in constant_def
    name = constant_def[KEY_NAME]
    t = constant_def[KEY_TYPE]
    value = constant_def[KEY_VALUE]
    
    q = '"' if t == YAML_TYPE_STR else ''
    
    f.write('const ' + yaml_type_to_cpp_type(t) + ' ' + name + ' = ' + q + str(value) + q + ';\n')
        
    
def generate_constants_header(constants_items, input_file_name):
    with open(os.path.join(TARGET_DIRECTORY, GEN_CONSTANTS_H), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f, input_file_name)
        
        guard = 'API_GEN_CONSTANTS';
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        
        f.write('#include <string>\n')
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        for constant_def in constants_items:
            write_constant_def(f, constant_def)

        f.write('\n' + DECL_DEFINE_PYBINDIND_CONSTANTS + ';\n\n')
        
        f.write(NAMESPACES_END + '\n\n')
        
        f.write('#endif // ' + guard + '\n\n')


def write_constant_binding(f, constant_def):
    assert KEY_NAME in constant_def
    assert KEY_TYPE in constant_def
    assert KEY_VALUE in constant_def
    name = constant_def[KEY_NAME]
    t = constant_def[KEY_TYPE]
    value = constant_def[KEY_VALUE]
    
    q = '"' if t == YAML_TYPE_STR else ''
    
    f.write('  m.attr("' + name + '") = py::' + yaml_type_to_pybind_type(t) + '(' + q + str(value) + q + ');\n')
    
        
def generate_constants_implementation(constants_items, input_file_name):
    with open(os.path.join(TARGET_DIRECTORY, GEN_CONSTANTS_CPP), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f, input_file_name)
        f.write(INCLUDE_API_COMMON_H + '\n')
        f.write('\n' + NAMESPACES_BEGIN + '\n\n')
        
        f.write(DECL_DEFINE_PYBINDIND_CONSTANTS + ' {\n')
        
        for constant_def in constants_items:
            write_constant_binding(f, constant_def)
        
        f.write('}\n\n')
        
        f.write(NAMESPACES_END + '\n\n')    


def generate_constants(constants_items, input_file_name):
    generate_constants_header(constants_items, input_file_name)
    generate_constants_implementation(constants_items, input_file_name)


def generate_data_classes(data_classes, input_file_name):
    assert type(data_classes) == dict
    for key, value in data_classes.items():
        if key != KEY_CONSTANTS:
            generate_class_files(data_classes, key, value, input_file_name)
        else:
            generate_constants(value, input_file_name)
    
        
def load_and_generate_data_classes():
    with open(DATA_CLASSES_FILE) as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        data_classes = yaml.load(file, Loader=yaml.FullLoader)
        if data_classes and type(data_classes) == dict:
            #print(data_classes)
            generate_data_classes(data_classes, DATA_CLASSES_FILE)
        else:
            print("Error while reading " + DATA_CLASSES_FILE + ".")
            sys.exit(1)

if __name__ == '__main__':
    load_and_generate_data_classes()
    
    