/*
 * mdl2datamodel.cpp
 *
 *  Created on: Mar 16, 2020
 *      Author: ahusar
 */
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "mdl2datamodel.h"
#include "datamodel_defines.h"
#include "json/json.h"

using namespace std;
using namespace Json;


// ------------- utilities -----------------------
static const char* get_sym_name(const sym_entry *s) {
  if (s == nullptr) {
    return "NULL SYMBOL";
  }
  else {
    return s->name;
  }
}

static bool has_flag(unsigned int flags, unsigned int one_flag) {
  assert(__builtin_popcount(one_flag) == 1); // check that we are testing just a single bit
  return (flags & one_flag) != 0;
}


// ------------- conversion methods -----------------------

bool Mdl2DataModel::convert_molecule(Value& molecule_list, species* spec) {

  const char* name = get_sym_name(spec->sym);
  if (strcmp(name, ALL_MOLECULES) == 0 ||
      strcmp(name, ALL_VOLUME_MOLECULES) == 0 ||
      strcmp(name, ALL_SURFACE_MOLECULES) == 0) {
    return true;
  }

  Value m;
  m[KEY_MOL_NAME] = name;
  m[KEY_DIFFUSION_CONSTANT] = spec->D;
  m[KEY_MOL_TYPE] = (has_flag(spec->flags, ON_GRID)) ? VALUE_MOL_TYPE_2D : VALUE_MOL_TYPE_3D;

  molecule_list.append(m);

  return true;
}


bool Mdl2DataModel::convert_molecule_list(Value& mcell) {

  Value& define_molecules = mcell[KEY_DEFINE_MOLECULES];
  Value& molecule_list = define_molecules[KEY_MOLECULE_LIST];

  for (int i = 0; i < s->n_species; i++) {
    species* spec = s->species_list[i];

    convert_molecule(molecule_list, spec);
  }

  DMUtil::json_add_version(define_molecules, JSON_DM_VERSION_1638);

  return true;
}


bool Mdl2DataModel::convert(const volume* state_, const char* output_file_name) {

  s = state_;

  cout << "Convert called\n";

  Value root;
  Value& mcell = root[KEY_MCELL];

  convert_molecule_list(mcell);

  Json::StreamWriterBuilder wbuilder;
  wbuilder["indentation"] = " ";
  std::string document = Json::writeString(wbuilder, root);

  // write result into a file
  ofstream myfile(output_file_name);
  if (myfile.is_open())
  {
    myfile << document;
    myfile.close();
  }
  else {
    cout << "Unable to open file " << output_file_name << " for writing.\n";
  }

  return true;
}


// C-language style entry point
bool convert_to_datamodel(const volume* state, const char* output_file_name) {
  Mdl2DataModel converter;
  return converter.convert(state, output_file_name);
}
