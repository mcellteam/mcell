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

enum class ReleaseShape {
  UNDEFINED = -1,  /* Not specified */
  SPHERICAL,       /* Volume enclosed by a sphere */
  // SHAPE_CUBIC,           /* Volume enclosed by a cube */ (might be supported, needs to be tested)
  // SHAPE_ELLIPTIC,        /* Volume enclosed by an ellipsoid */ (might be supported, needs to be tested)
  // SHAPE_RECTANGULAR,     /* Volume enclosed by a rect. solid */ (might be supported, needs to be tested)
  SPHERICAL_SHELL, /* Surface of a sphere */ // not tested yet
  REGION,          /* Inside/on the surface of an arbitrary region */
  // SHAPE_LIST             /* Individiaul mol. placement by list */
};

/**
 * Release molecules according to the settings.
 */
class ReleaseEvent: public BaseEvent {
public:
  ReleaseEvent(World* world_) :
    BaseEvent(EVENT_TYPE_INDEX_RELEASE),
    release_site_name(NAME_INVALID),
    species_id(SPECIES_ID_INVALID),
    release_number(0),
    orientation(ORIENTATION_NONE),
    release_shape(ReleaseShape::UNDEFINED),
    world(world_) {
  }
  virtual ~ReleaseEvent() {}

  virtual void step();
  virtual void dump(const std::string indent);

public:
  std::string release_site_name; // name of releaser site from which was this event created

  species_id_t species_id;
  uint release_number; // number of molecules to release

  orientation_t orientation;

  ReleaseShape release_shape; /* Release Shape Flags: controls shape over which to release (enum release_shape_t) */

  // SHAPE_SPHERICAL - only volume molecules
  vec3_t location;
  vec3_t diameter; /* x,y,z diameter for geometrical release shapes */

  // SHAPE_REGION
  // for surface molecule releases
  std::vector<CummAreaPWallIndexPair> cum_area_and_pwall_index_pairs;

  // for volume molecule releases into a region
  std::string region_name; // name of the region into which we should release the
  vec3_t region_llf; // note: this is fully specified by the region above, maybe remove in the future
  vec3_t region_urb; // note: this is fully specified by the region above as well


  std::string release_pattern_name; // unused

private:
  World* world;

private:
  uint calculate_number_to_release();

  // for surface molecule releases
  void place_single_molecule_onto_grid(Partition& p, Wall& wall, uint tile_index);
  void release_onto_regions(uint computed_release_number);

  // for volume molecule releases into a region
  void release_inside_regions(uint computed_release_number);

  // for volume molecule releases
  void release_ellipsoid_or_rectcuboid(uint computed_release_number);

};

} // namespace mcell


#endif // SRC4_RELEASE_EVENT_H_
