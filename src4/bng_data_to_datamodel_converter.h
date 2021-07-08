/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_BNGDATA_TO_DATAMODEL_CONVERTER_H_
#define SRC4_BNGDATA_TO_DATAMODEL_CONVERTER_H_

#include <cassert>
#include "defines.h"

namespace Json {
class Value;
}

namespace BNG {
class BNGEngine;
class RxnRule;
class ElemMolType;
}

namespace MCell {

class World;

/**
 * We want to keep the BNG independent,
 * so this is an extra too to convert the GND data to the
 * cellblender datamodel.
 */
class BngDataToDatamodelConverter {
public:
  BngDataToDatamodelConverter();

  // does nothing for now, there will be changes in BNG data and
  // converting species is not needed at this point
  void to_data_model(const World* world_, Json::Value& mcell_node, const bool only_for_viz);

private:
  void reset();

  Vec3 get_next_color();
  void convert_molecules(Json::Value& mcell_node);
  void convert_single_mol_type(const BNG::ElemMolType& s, Json::Value& molecule_node);

  void convert_single_rxn_rule(const BNG::RxnRule& r, Json::Value& species_node);
  std::string get_surface_class_name(const BNG::RxnRule& r);
  void convert_single_surface_class(const BNG::RxnRule& r, Json::Value& rxn_node);
  void convert_rxns(Json::Value& mcell_node);

  const World* world;
  const BNG::BNGEngine* bng_engine;

  uint next_color_index;
  uint rxn_counter;
  std::set<std::string> processed_surface_classes;
  std::vector<Vec3> colors;

  bool conversion_failed;
};

} // namespace MCell

#endif // SRC4_BNGDATA_TO_DATAMODEL_CONVERTER_H_
