/******************************************************************************
 *
 * Copyright (C) 2019 by
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

#ifndef SRC4_RELEASE_EVENT_H_
#define SRC4_RELEASE_EVENT_H_

#include <vector>

#include "base_event.h"

namespace MCell {

class Partition;
class Wall;
class Grid;


/**
 * Release molecules according to the settings.
 */
class ReleaseEvent: public BaseEvent {
public:
  ReleaseEvent(World* world_) :
    BaseEvent(EVENT_TYPE_INDEX_RELEASE),
    species_id(SPECIES_ID_INVALID),
    release_number(0),
    release_site_name(NAME_INVALID),
    orientation(ORIENTATION_NONE),
    release_shape(SHAPE_SPHERICAL),
    world(world_) {
  }
  virtual ~ReleaseEvent() {}

  virtual void step();
  virtual void dump(const std::string indent);

public:
  vec3_t location;
  species_id_t species_id;
  uint release_number; // number of molecules to release

  std::string release_pattern_name;

  // ---------------------------------- release location ----------------------------------

  std::string release_site_name;
  orientation_t orientation;
  int8_t release_shape; /* Release Shape Flags: controls shape over which to release (enum release_shape_t) */
  vec3_t diameter; /* x,y,z diameter for geometrical release shapes */

  // ---------------------------------- surface release info ----------------------------------

  // for surface molecule releases
  std::vector<CummAreaPWallIndexPair> cum_area_and_pwall_index_pairs;

  // do I need objects?
  // walls?
  std::vector<wall_index_t> wall_indices_t; // need to find the correct partition?


  World* world;

private:
  uint calculate_number_to_release();

  void place_single_molecule_onto_grid(Partition& p, Wall& wall, uint tile_index);
  void release_onto_regions(uint computed_release_number);

  void release_ellipsoid_or_rectcuboid(uint computed_release_number);

};

} // namespace mcell


#endif // SRC4_RELEASE_EVENT_H_