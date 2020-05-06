/******************************************************************************
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
class Species;
class RxnRule;
}

namespace MCell {

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
  void to_data_model(const BNG::BNGEngine& bng_engine, Json::Value& mcell_node);

private:
  void reset();

  Vec3 get_next_color();
  void convert_species(Json::Value& mcell_node);
  void convert_single_species(const BNG::Species& s, Json::Value& species_node);

  void convert_single_rxn_rule(const BNG::RxnRule& r, Json::Value& species_node);
  void convert_rxns(Json::Value& mcell_node);

  const BNG::BNGEngine* bng_engine;

  uint next_color_index;
  std::vector<Vec3> colors;

  bool conversion_failed;
};

} // namespace MCell

#endif // SRC4_BNGDATA_TO_DATAMODEL_CONVERTER_H_
