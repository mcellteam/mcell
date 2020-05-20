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

#ifndef LIBMCELL_API_MCELL4_CONVERTER_H_
#define LIBMCELL_API_MCELL4_CONVERTER_H_

#include <vector>
#include "defines.h"

namespace BNG {
class ComponentType;
class ComponentInstance;
class MolType;
class MolInstance;
class CplxInstance;
}

namespace MCell {

class World;
class Partition;
class GeometryObject;

namespace API {


class Model;
class ComponentType;
class ComponentInstance;
class MoleculeType;
class MoleculeInstance;
class ComplexInstance;
class GeometryObject;

class MCell4Converter {
public:
  // throws exception if anything went wrong
  // modifies model as well where it stores information for cases
  // when a value such a reaction rate was updated by the user
  void convert(Model* model_, World* world_);

private:

  void convert_simulation_setup();

  void convert_species();

  BNG::ComponentType convert_component_type(API::ComponentType& ct);
  BNG::ComponentInstance convert_component_instance(API::ComponentInstance& ci);
  BNG::MolType convert_molecule_type(API::MoleculeType& mt);
  BNG::MolInstance convert_molecule_instance(API::MoleculeInstance& mi);
  BNG::CplxInstance convert_complex_instance(API::ComplexInstance& inst);
  void convert_rxns();

  wall_index_t convert_wall_and_add_to_geom_object(
      const API::GeometryObject& src_obj, const uint side,
      MCell::Partition& p, MCell::GeometryObject& dst_obj);

  void convert_geometry_objects();

  void convert_release_events();

  void convert_viz_output_events();

  Model* model;
  World* world;
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_MCELL4_CONVERTER_H_ */
