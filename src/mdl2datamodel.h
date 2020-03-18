/*
 * mdl2datamodel.h
 *
 *  Created on: Mar 16, 2020
 *      Author: ahusar
 */

#ifndef SRC_MDL2DATAMODEL_H_
#define SRC_MDL2DATAMODEL_H_

#include "mcell_structs.h"
#include "json/json.h"

class Mdl2DataModel {
public:
  bool convert(const volume* state_, const char* output_file_name);

private:
  void add_version(Json::Value& define_molecules, const char* ver);
  bool convert_molecule(Json::Value& molecule_list, species* spec);
  bool convert_molecule_list(Json::Value& mcell);

  // state initialized when convert
  const volume* s;
};

bool convert_to_datamodel(const volume* state, const char* output_file_name);

#endif /* SRC_MDL2DATAMODEL_H_ */
