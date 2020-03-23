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

# TODO: change 'items' to 'attributes'

import sys
import os
import yaml
import re
from datetime import datetime

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

KEY_SUPERCLASS = 'superclass'

KEY_ITEMS = 'items'
KEY_NAME = 'name'
ATTR_NAME_NAME = 'name' # attrribute with name 'name' is already defined in BaseDataClass
KEY_TYPE = 'type'
KEY_DEFAULT = 'default'

KEY_METHODS = 'methods'
KEY_PARAMS = 'params'
KEY_RETURN_TYPE = 'return_type'

YAML_TYPE_FLOAT = 'float'
YAML_TYPE_STR = 'str'
YAML_TYPE_INT = 'int'
YAML_TYPE_VEC2 = 'Vec2'
YAML_TYPE_VEC3 = 'Vec3'

CPP_TYPE_FLOAT = 'float_t'
CPP_TYPE_STR = 'std::string'
CPP_TYPE_INT = 'int'
CPP_TYPE_VEC2 = 'Vec2'
CPP_TYPE_VEC3 = 'Vec3'

CPP_REFERENCE_TYPES = [CPP_TYPE_STR, CPP_TYPE_VEC2, CPP_TYPE_VEC3]

UNSET_VALUE_FLOAT = 'FLT_UNSET'
UNSET_VALUE_STR = 'STR_UNSET'
UNSET_VALUE_INT = 'INT_UNSET'
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
RET_TYPE_TO_STR = 'std::string'
DECL_TO_STR = 'to_str() const'
KEYWORD_OVERRIDE = 'override'
  
CLASS_NAME_ATTR = 'class_name'

INCLUDE_API_MCELL_H = '#include "../api/mcell.h"'
NAMESPACES_BEGIN = 'namespace MCell {\nnamespace API {'
NAMESPACES_END = '} // namespace API\n} // namespace MCell'


def get_underscored(class_name):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', class_name).lower()

def get_header_guard_name(class_name):
    return GUARD_PREFIX + get_underscored(class_name).upper() + GUARD_SUFFIX

def get_class_file_name(class_name, extension):
    return os.path.join(TARGET_DIRECTORY, GEN_PREFIX + get_underscored(class_name) + '.' + extension)

def yaml_type_to_cpp_type(t):
    assert len(t) >= 1
    if t == YAML_TYPE_FLOAT:
        return CPP_TYPE_FLOAT
    elif t == YAML_TYPE_STR:
        return CPP_TYPE_STR
    elif t == YAML_TYPE_INT:
        return CPP_TYPE_INT
    elif t == YAML_TYPE_VEC2:
        return CPP_TYPE_VEC2
    elif t == YAML_TYPE_VEC3:
        return CPP_TYPE_VEC3
    else:
        return t + '*' # everything else is passed as a pointer
    

def get_type_as_ref_param(attr):
    assert KEY_TYPE in attr
    
    yaml_type = attr[KEY_TYPE]
    cpp_type = yaml_type_to_cpp_type(yaml_type)
    
    res = cpp_type
    if cpp_type in CPP_REFERENCE_TYPES:
        res += '&'
    
    return res


def get_unset_value(attr):
    assert KEY_TYPE in attr
    t = attr[KEY_TYPE]
    if t == YAML_TYPE_FLOAT:
        return UNSET_VALUE_FLOAT
    elif t == YAML_TYPE_STR:
        return UNSET_VALUE_STR
    elif t == YAML_TYPE_INT:
        return UNSET_VALUE_INT
    elif t == YAML_TYPE_VEC2:
        return UNSET_VALUE_VEC2
    elif t == YAML_TYPE_VEC3:
        return UNSET_VALUE_VEC3
    else:
        return UNSET_VALUE_PTR


def is_cpp_ptr_type(t):
    assert len(t) >= 1
    return t[-1] == '*'


def is_yaml_ptr_type(t):
    return is_cpp_ptr_type(yaml_type_to_cpp_type(t))

    
def is_cpp_ptr_or_ref_type(t):
    return is_cpp_ptr_type(t) or t in CPP_REFERENCE_TYPES


def write_generated_notice(f, input_file_name):
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M")
    f.write('// This file was generated automatically on ' + date_time + ' from ' + '\'' + input_file_name + '\'\n\n')


def write_ctor_define(f, class_name, items):
    f.write('#define ' + get_underscored(class_name).upper() + CTOR_SUFFIX + '() \\\n')
    f.write('    ' + class_name + '( \\\n')

    # ctor parameters    
    num_items = len(items)
    for i in range(num_items):
        attr = items[i]
        assert KEY_NAME in attr
        name = attr[KEY_NAME]

        f.write('        const ' + get_type_as_ref_param(attr) + ' ' + name + '_')
        
        if KEY_DEFAULT in attr:
            f.write(' = ' + get_unset_value(attr))
        
        if i != num_items - 1:
            f.write(',')
        f.write(' \\\n')

    f.write('    ) { \\\n')
    
    # initialization code
    f.write('      ' + CLASS_NAME_ATTR + ' = "' + class_name + '"; \\\n')
    for i in range(num_items):
        assert KEY_NAME in items[i] 
        attr_name = items[i][KEY_NAME]
        f.write('      ' + attr_name + ' = ' + attr_name + '_; \\\n')
    f.write('    }\n\n')    
    

def write_attr_with_get_set(f, attr):
    assert KEY_NAME in attr
    assert KEY_TYPE in attr
    
    name = attr[KEY_NAME]
    
    # skip attribute 'name'
    if name == ATTR_NAME_NAME:
        return False
    
    yaml_type = attr[KEY_TYPE]
    cpp_type = yaml_type_to_cpp_type(yaml_type)
    
    decl_const = 'const ' if is_cpp_ptr_type(cpp_type) else ''
    decl_type = decl_const + cpp_type
    
    # decl
    f.write('  ' + decl_type + ' ' + name + ';\n')
    
    # setter
    f.write('  virtual void set_' + name + '(const ' + get_type_as_ref_param(attr) + ' new_' + name + '_) {\n')
    f.write('    ' + name + ' = new_' + name + '_;\n')
    f.write('  }\n') 

    # getter
    ret_type_const = 'const ' if is_cpp_ptr_or_ref_type(cpp_type) else ''
    
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
            f.write('const ' + t + ' ' + p[KEY_NAME])
            if i != params_cnt - 1:
                f.write(', ')
                
    f.write(')')
    
        
def write_method_declaration(f, method):
    f.write('  virtual ')
    write_method_signature(f, method)
    f.write(' = 0;\n')
    

def based_on_base_data_class(class_def):
    if KEY_SUPERCLASS in class_def:
        assert class_def[KEY_SUPERCLASS] == BASE_DATA_CLASS # the only supported superclass now
        return True
    else:
        return False
    
    
def write_gen_class(f, class_name, class_def):
    f.write('class ' + GEN_CLASS_PREFIX + class_name)
    
    if based_on_base_data_class(class_def):
        f.write(': public ' + BASE_DATA_CLASS)
        
    f.write( ' {\n')
    f.write('public:\n')
    
    if based_on_base_data_class(class_def):
        f.write('  ' + RET_TYPE_CHECK_SEMANTICS + ' ' + DECL_CHECK_SEMANTICS + ' ' + KEYWORD_OVERRIDE + ';\n')
        f.write('  ' + RET_TYPE_TO_STR + ' ' +DECL_TO_STR + ' ' + KEYWORD_OVERRIDE + ';\n\n')
        
    f.write('  // --- attributes ---\n')
    items = class_def[KEY_ITEMS]
    for attr in items:
        written = write_attr_with_get_set(f, attr)
        if written :
            f.write('\n')

    f.write('  // --- methods ---\n')
    methods = class_def[KEY_METHODS]
    for m in methods:
        write_method_declaration(f, m)
        
    f.write('}; // ' + GEN_CLASS_PREFIX + class_name + '\n\n')


def write_define_binding_decl(f, class_name):
    f.write('void define_pybinding_' + class_name + '(py::module& m)')           


def generate_class_header(class_name, class_def, input_file_name):
    with open(get_class_file_name(class_name, EXT_H), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f, input_file_name)
        
        guard = get_header_guard_name(class_name);
        f.write('#ifndef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        f.write(INCLUDE_API_MCELL_H + '\n\n')
        f.write(NAMESPACES_BEGIN + '\n\n')
        
        if based_on_base_data_class(class_def):
            write_ctor_define(f, class_name, class_def[KEY_ITEMS])
        
        write_gen_class(f, class_name, class_def)
        
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
    f.write(RET_TYPE_CHECK_SEMANTICS + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_CHECK_SEMANTICS + '{\n') 
    for attr in items:
        if KEY_DEFAULT not in attr:
            write_is_set_check(f, attr[KEY_NAME])
    
    f.write('  return SemRes::OK;\n')
    f.write('}\n\n')    
    
        
def write_to_str_implemetation(f, class_name, items):
    f.write(RET_TYPE_TO_STR + ' ' + GEN_CLASS_PREFIX + class_name + '::' + DECL_TO_STR + '{\n')
    f.write('  std::stringstream ss;\n')
    f.write('  ss << get_object_name() << ": " <<\n')
    
    num_attrs = len(items) 
    for i in range(num_attrs):
        name = items[i][KEY_NAME]
        f.write('      "' + name + '=" << ')
        
        type = items[i][KEY_TYPE]
        if is_yaml_ptr_type(type):
            f.write('((' + name + ' != nullptr) ? ' + name + '->to_str() : "null" )')
        else:
            f.write(name)

        if i != num_attrs - 1:
            f.write(' << ", " <<\n')
            
    f.write(';\n')
    
    f.write('  return ss.str();\n')
    f.write('}\n\n')                


def write_pybind11_bindings(f, class_name, class_def):
    items = class_def[KEY_ITEMS]
    
    write_define_binding_decl(f, class_name)
    f.write(' {\n')
    f.write('  py::class_<' + class_name + '>(m, "' + class_name + '")\n')
    f.write('      .def(\n')
    f.write('          py::init<\n')

    # init operands
    num_items = len(items)
    for i in range(num_items):
        attr = items[i]
        f.write('            const ' + get_type_as_ref_param(attr))
        if i != num_items - 1:
            f.write(',\n')
    f.write('\n')
    f.write('          >(),\n')
        
    # init argument names and default values
    for i in range(num_items):
        attr = items[i]
        name = attr[KEY_NAME]
        f.write('          py::arg("' + name + '")')
        if KEY_DEFAULT in attr:
            f.write(' = ' + get_unset_value(attr))
        if i != num_items - 1:
            f.write(',\n')
    f.write('\n')            
    f.write('        )\n')            
    
    # common methods
    if based_on_base_data_class(class_def):
        f.write('      .def("check_semantics", &' + class_name + '::check_semantics_cerr)\n')
        f.write('      .def("__str__", &' + class_name + '::to_str)\n')
        
    # dump is required to be implemented, TODO: to_str as well?
    f.write('      .def("dump", &' + class_name + '::dump)\n')
    
    # properties
    for i in range(num_items):
        name = items[i][KEY_NAME]
        f.write('      .def_property("' + name + '", &' + class_name + '::get_' + name + ', &' + class_name + '::set_' + name + ')\n')
    f.write('    ;\n')
    f.write('}\n\n')
                
            
def generate_class_implementation_and_bindings(class_name, class_def, input_file_name):
    with open(get_class_file_name(class_name, EXT_CPP), 'w') as f:
        f.write(COPYRIGHT)
        write_generated_notice(f, input_file_name)
        
        f.write('#include <sstream>\n')
        f.write(INCLUDE_API_MCELL_H + '\n')
        
        f.write(NAMESPACES_BEGIN + '\n\n')
        
        if based_on_base_data_class(class_def):
            items = class_def[KEY_ITEMS]
            write_check_semantics_implemetation(f, class_name, items)
            write_to_str_implemetation(f, class_name, items)
        
        write_pybind11_bindings(f, class_name, class_def)
        
        f.write(NAMESPACES_END + '\n\n')
        


def generate_class_files(class_name, class_def, input_file_name):
    # we need items and methods to be present 
    if KEY_ITEMS not in class_def:
        class_def[KEY_ITEMS] = []
        
    if KEY_METHODS not in class_def:
        class_def[KEY_METHODS] = []
    
    generate_class_header(class_name, class_def, input_file_name)
    generate_class_implementation_and_bindings(class_name, class_def, input_file_name)
    

def generate_data_classes(data_classes, input_file_name):
    assert type(data_classes) == dict
    for key, value in data_classes.items():
        generate_class_files(key, value, input_file_name)
    
    
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
    
    