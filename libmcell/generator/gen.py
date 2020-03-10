#!/usr/bin/env python3

import sys
import os
import yaml
import re

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
******************************************************************************/
\n
"""

KEY_ITEMS = 'items'
KEY_NAME = 'name'
KEY_TYPE = 'type'
KEY_DEFAULT = 'default'

YAML_TYPE_FLOAT = 'float'
YAML_TYPE_STR = 'str'
YAML_TYPE_INT = 'int'
YAML_TYPE_VEC2 = 'Vec2'
YAML_TYPE_VEC3 = 'Vec3'
YAML_TYPE_PTR = 'ptr' # not specific, original type name must be used

CPP_TYPE_FLOAT = 'float_t'
CPP_TYPE_STR = 'std::string'
CPP_TYPE_INT = 'int'
CPP_TYPE_VEC2 = 'Vec2'
CPP_TYPE_VEC3 = 'Vec3'
CPP_TYPE_PTR = 'ptr' # not specific, original type name must be used

CPP_REFERENCE_TYPES = [CPP_TYPE_STR, CPP_TYPE_VEC2, CPP_TYPE_VEC3]

UNSET_VALUE_FLOAT = 'FLT_UNSET'
UNSET_VALUE_STR = 'STR_UNSET'
UNSET_VALUE_INT = 'INT_UNSET'
UNSET_VALUE_VEC2 = 'VEC2_UNSET'
UNSET_VALUE_VEC3 = 'VEC3_UNSET'
UNSET_VALUE_PTR = 'nullptr'

GEN_PREFIX = 'gen_'
GUARD_PREFIX = 'API_GEN_'
GUARD_SUFFIX = 'API_GEN_'
CTOR_SUFFIX = '_CTOR'
EXT_CPP = 'cpp'
EXT_H = 'h'

CLASS_NAME_ATTR = 'class_name'

INCLUDE_API_MCELL_H = '#include "../api/mcell.h'
NAMESPACES_BEGIN = 'namespace MCell {\nnamespace API {'
NAMESPACES_END = '} // namespace API\n} // namespace MCell'


def get_underscored(class_name):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', class_name).lower()

def get_header_guard_name(class_name):
    return GUARD_PREFIX + get_underscored(class_name).upper() + GUARD_SUFFIX

def get_class_file_name(class_name, extension):
    return GEN_PREFIX + get_underscored(class_name) + '.' + extension

def yam_type_to_cpp_type(t):
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
    elif t[-1] == '*':
        return CPP_TYPE_PTR
    else:
        print("Unexpected type '" + t + "'")
        assert False
    

def get_as_ref_param(attr):
    assert KEY_NAME in attr
    assert KEY_TYPE in attr
    
    yaml_type = attr[KEY_TYPE]
    cpp_type = yam_type_to_cpp_type(yaml_type)
    
    res = ''
    if cpp_type != CPP_TYPE_PTR:
        res += cpp_type
    else:
        res += yaml_type
    
    if t in CPP_REFERENCE_TYPES:
        res += '&'
    
    res += ' ' + attr[KEY_NAME]
    return res


def get_unset_value(attr):
    assert KEY_TYPE in attr
    t = attr[KEY_TYPE]
    assert len(t) >= 1
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
    elif t[-1] == '*':
        return UNSET_VALUE_PTR
    else:
        print("Unexpected type '" + t + "'")
        assert False


def write_ctor_define(f, class_name, items):
    f.write('#define ' + get_underscored().upper() + CTOR_SUFFIX + '() \\')
    f.write('\t' + class_name + '( \\')

    # ctor parameters    
    num_items = len(items)
    for i in range(num_items):
        attr = items[i]
        f.write('\t\tconst ' + get_as_ref_param(attr) + '_')
        
        if KEY_DEFAULT in attr:
            f.write(' = ' + get_unset_value(attr))
        
        if i != num_items - 1:
            f.write(',')
        f.write(' \\')

    f.write('\t) { \\')
    
    # initialization code
    f.write('\t' + CLASS_NAME_ATTR + ' = "' + class_name + '"; \\')
    for i in range(num_items):
        assert KEY_NAME in items[i] 
        attr_name = items[i][KEY_NAME]
        f.write('\t  ' + attr_name + ' = ' + attr_name + '_; \\')
    f.write('\t}\n')    
    

def generate_class_header(class_name, items):
    with open(get_class_file_name(class_name, EXT_H), 'w') as f:
        f.write(COPYRIGHT)
        guard = get_header_guard_name();
        f.write('#ifdef ' + guard + '\n')
        f.write('#define ' + guard + '\n\n')
        f.write(INCLUDE_API_MCELL_H + '\n\n')
        f.write(NAMESPACES_BEGIN)
        
        write_ctor_define(f, class_name, items)
        
        f.write(NAMESPACES_END)
        f.write('#endif // ' + guard + '\n')
    

def generate_class_files(class_name, contents):
    assert KEY_ITEMS in contents
    items = contents[KEY_ITEMS]
    
    generate_class_header(class_name, items)
    

def generate_data_classes(data_classes):
    assert type(data_classes) == dict
    for key, value in data_classes.items():
        generate_class_files(key, value)
    
    
def load_and_generate_data_classes():
    with open(DATA_CLASSES_FILE) as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        data_classes = yaml.load(file, Loader=yaml.FullLoader)
        if data_classes and type(data_classes) == dict:
            print(data_classes)
            generate_data_classes(data_classes)
        else:
            print("Error while reading " + DATA_CLASSES_FILE + ".")
            sys.exit(1)

if __name__ == '__main__':
    load_and_generate_data_classes()
    
    