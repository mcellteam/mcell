"""
Copyright (C) 2021 by
The Salk Institute for Biological Studies 

Use of this source code is governed by an MIT-style
license that can be found in the LICENSE file or at
https://opensource.org/licenses/MIT.
"""

import sys
import os
import yaml

from constants import *
from gen import indent_and_fix_rst_chars, yaml_type_to_py_type, get_default_or_unset_value_py


def cat_to_title(cat):
    if cat == CATEGORY_CONSTANTS:
        return 'Enums and Constants'
    else:
        return cat.replace('_', ' ').capitalize()
        
def write_cat_label(f, cat):
    f.write('.. _api-' + cat + ':\n\n')
            
        
def gen_example_links(base_links):
    split_links = base_links.strip().split()
    n = len(split_links)
    if n == 0:
        return ''        
    res = 'Example' + ('' if n == 1 else 's') + ': '
    
    for l in split_links:
        name = os.path.basename(os.path.dirname(l)) + '/' + os.path.basename(l)
        res += '`' + name + ' <' + EXAMPLES_BASE_URL + l + '>`_ '
    
    return res 


def write_h4(f, text, name, class_name):
    
    f.write('.. _' + class_name + '__' + name + ':\n\n')
    
    f.write(text + '\n')
    f.write('-' * len(text) + '\n\n')


def get_method_declaration(method):
    res = method[KEY_NAME] + ' ('
    
    if KEY_PARAMS in method:
        num_params = len(method[KEY_PARAMS])
        for i in range(num_params):
            param = method[KEY_PARAMS][i]
            t = yaml_type_to_py_type(param[KEY_TYPE])
            res += param[KEY_NAME] + ': ' + t
            
            if KEY_DEFAULT in param:
                res += '=' + get_default_or_unset_value_py(param)
                
            if i != num_params - 1:
                res += ', '

    res += ')'

    if KEY_RETURN_TYPE in method:
        res += ' -> ' + yaml_type_to_py_type(method[KEY_RETURN_TYPE])
        
    return res
            

def generate_class_documentation(f, class_name, class_def):
    f.write(class_name + '\n' + '='*len(class_name) + '\n\n')
    
    if KEY_DOC in class_def:
        f.write(class_def[KEY_DOC].strip() + '\n\n')
        
    if KEY_EXAMPLES in class_def:
        f.write(gen_example_links(class_def[KEY_EXAMPLES]) + '\n\n')
        
        
    if KEY_ITEMS in class_def and class_def[KEY_ITEMS]:
        f.write('Attributes:\n' + '*'*len('Attributes:') + '\n')
        num_items = len(class_def[KEY_ITEMS])
        for item in class_def[KEY_ITEMS]:
            t = yaml_type_to_py_type(item[KEY_TYPE])
            
            header = item[KEY_NAME] + ': ' + t
            write_h4(f, header, item[KEY_NAME], class_name)
            
            if KEY_DOC in item and item[KEY_DOC]:
                f.write('  | ' + indent_and_fix_rst_chars(item[KEY_DOC].strip(), '  | ') + '\n')

            if KEY_DEFAULT in item:
                f.write('  | - default argument value in constructor: ' + get_default_or_unset_value_py(item))
            f.write('\n')
                
            if KEY_EXAMPLES in item:
                f.write('\n  | ' + gen_example_links(item[KEY_EXAMPLES]) + '\n\n')
                
            f.write('\n')
        
    if KEY_METHODS in class_def and class_def[KEY_METHODS]:
        f.write('\nMethods:\n' + '*'*len('nMethods:') + '\n')
        for method in class_def[KEY_METHODS]:
            method_name = method[KEY_NAME]
            
            header = get_method_declaration(method)
            write_h4(f, header, method_name, class_name)
            
            if KEY_DOC in method:
                f.write('\n  | ' + indent_and_fix_rst_chars(method[KEY_DOC].strip(), '  | ') + '\n\n')

            if KEY_PARAMS in method:
                num_params = len(method[KEY_PARAMS])
                for param in method[KEY_PARAMS]:
                    t = yaml_type_to_py_type(param[KEY_TYPE])
                    f.write('* | ' + param[KEY_NAME] + ': ' + t)
                    if KEY_DEFAULT in param:
                        f.write(' = ' + get_default_or_unset_value_py(param))
                    if KEY_DOC in param:
                        f.write('\n  | ' + indent_and_fix_rst_chars(param[KEY_DOC].strip(), '  | ') + '\n\n')
                    else:
                        f.write('\n')

            if KEY_EXAMPLES in method:
                f.write('  | ' + gen_example_links(method[KEY_EXAMPLES]) + '\n\n')
                
            f.write('\n')
        f.write('\n')
    
            
def generate_documentation(data_classes):
    
    # generate constants
    with open(os.path.join(DOC_DIRECTORY, CATEGORY_CONSTANTS + EXT_RST), 'w') as f:
        write_cat_label(f, CATEGORY_CONSTANTS)
        f.write(
            '*******************\n' +
            cat_to_title(CATEGORY_CONSTANTS) + '\n' +
            '*******************\n\n'
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
        
        
    # then generate classes into files by category
    for cat in CATEGORIES:
        if cat == CATEGORY_CONSTANTS:
            continue
        input_file = cat + EXT_RST
        with open(os.path.join(DOC_DIRECTORY, input_file), 'w') as f:
            write_cat_label(f, cat)
            cat_name = cat_to_title(cat)
            f.write('*'*len(cat_name) + '\n' + cat_name + '\n' + '*'*len(cat_name) + '\n')
                 
            for key, value in sorted(data_classes.items()):
                if key != KEY_CONSTANTS and key != KEY_ENUMS and value[KEY_CATEGORY] == cat:
                    generate_class_documentation(f, key, value)

    # and generate api.rst file
    with open(os.path.join(DOC_DIRECTORY, API_RST), 'w') as f:

        title = 'Python API Reference'
        f.write(
            title + '\n' +
            '='*len(title) + '\n\n'
        )
        
        f.write(
            '.. toctree::\n'
            '    :maxdepth: 2\n'
            '    :hidden:\n'
            '    :caption: Contents\n\n'
        )
        
        for cat in CATEGORIES:
            f.write('    ' + cat + '\n')
        
        f.write('\nThis section contains automatically generated documentation on Python classes, enums, '
            'and constants provided by MCell.\n\n')
        
        for cat in CATEGORIES:
            f.write('- :ref:`api-' + cat + '`\n')
            