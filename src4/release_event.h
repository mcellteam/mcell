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

namespace Json {
class Value;
}

namespace MCell {

class Partition;
class Wall;
class Grid;

// TODO: replace with defines from API
enum class ReleaseShape {
  UNDEFINED = -1,  /* Not specified */
  SPHERICAL,       /* Volume enclosed by a sphere */
  // CUBIC,           /* Volume enclosed by a cube */ (might be supported, needs to be tested)
  // ELLIPTIC,        /* Volume enclosed by an ellipsoid */ (might be supported, needs to be tested)
  // RECTANGULAR,     /* Volume enclosed by a rect. solid */ (might be supported, needs to be tested)
  SPHERICAL_SHELL, /* Surface of a sphere */ // not tested yet
  REGION,          /* Inside/on the surface of an arbitrary region */
  // LIST             /* Individiaul mol. placement by list */
};

enum class ReleaseNumberMethod {
  Invalid,
  ConstNum,
  GaussNum,
  VolNum,
  ConcNum,
  DensityNum
};


// TODO: unify with API::RegionNodeType?
enum class RegionExprOperator {
  Invalid,
  Union,
  Intersect,
  Difference,
  Leaf
};

const int NUMBER_OF_TRAINS_UNLIMITED = -1;

class RegionExprNode {
public:
  RegionExprNode()
    : op(RegionExprOperator::Invalid), left(nullptr), right(nullptr) {
  }

  RegionExprOperator op;

  std::string region_name; // name of the region into which we should release the molecules

  RegionExprNode* left;
  RegionExprNode* right;

  void dump() const; // does not print any newlines
  std::string to_string(const bool for_datamodel = false) const;
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
    actual_release_time(TIME_INVALID),
    release_number_method(ReleaseNumberMethod::Invalid),
    release_number(UINT_INVALID),
    concentration(FLT_INVALID),
    orientation(ORIENTATION_NONE),
    release_shape(ReleaseShape::UNDEFINED),
    location(FLT_INVALID),
    diameter(FLT_INVALID),
    region_expr_root(nullptr),
    // the default values for release pattern are such that there is
    // just one release
    delay(0),
    number_of_trains(1),
    train_interval(EPS),
    train_duration(EPS),
    release_interval(FLT_GIGANTIC),
    current_train_from_0(0),
    current_release_in_train_from_0(0),
    world(world_)
     {
  }
  virtual ~ReleaseEvent();

  void step() override;

  // release events must be sorted by the actual release time as well
  bool needs_secondary_ordering() override {
    return true;
  }

  float_t get_secondary_ordering_value() override {
    assert(actual_release_time != TIME_INVALID);
    return actual_release_time;
  }

  // handles release patterns,
  // must be called before the first insertion into the scheduler,
  // when no release pattern is set, simply sets event_time to 0
  // and on the second call it returns false
  bool update_event_time_for_next_scheduled_time();

  void dump(const std::string indent) const override;
  void to_data_model(Json::Value& mcell_node) const;

public:
  std::string release_site_name; // name of releaser site from which was this event created

  species_id_t species_id;

  float_t actual_release_time;

  ReleaseNumberMethod release_number_method; // specifies what does the release_number mean
  uint release_number; // number of molecules to release
  float_t concentration;

  orientation_t orientation;

  ReleaseShape release_shape; /* Release Shape Flags: controls shape over which to release (enum release_shape_t) */

  // SHAPE_SPHERICAL - only volume molecules
  Vec3 location;
  Vec3 diameter; /* x,y,z diameter for geometrical release shapes */

  // SHAPE_REGION
  // for surface molecule releases

  // limited initialization for pymcell4
  // return true if initialization passed
  bool initialize_walls_for_release();

  // initialized from mcell3 state or in initialize_walls_for_release()
  // TODO: replace with some reasonale strcucture
  std::vector<CummAreaPWallIndexPair> cumm_area_and_pwall_index_pairs;


  // constructor and container for all region expr nodes
  RegionExprNode* create_new_region_expr_node_leaf(const std::string region_name);
  RegionExprNode* create_new_region_expr_node_op(const RegionExprOperator op, RegionExprNode* left, RegionExprNode* right);
  std::vector<RegionExprNode*> all_region_expr_nodes;

  // only when release_shape is SHAPE_REGION
  RegionExprNode* region_expr_root;

  // for volume molecule releases into a region

  Vec3 region_llf; // note: this is fully specified by the region above, maybe remove in the future
  Vec3 region_urb; // note: this is fully specified by the region above as well


  std::string release_pattern_name;

  // --- release pattern information ---
  float_t delay;
  int number_of_trains; // -1 means that the number of trains is unlimited
  float_t train_interval;
  float_t train_duration;
  float_t release_interval;


private:
  // both values are initialized to 0 and counted from 0
  int current_train_from_0;
  int current_release_in_train_from_0;

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

  float_t get_release_delay_time() const {
    if (cmp_eq(actual_release_time, event_time)) {
      return 0; // same as event time
    }
    else {
      return actual_release_time - event_time;
    }
  }

  int get_num_releases_per_train() const {
    assert(release_interval != 0);

    if (train_duration > EPS) {
      // -EPS is needed to deal with precision issues even when we
      // are strictly (<) comparing current and end time
      return ceil_f((train_duration - EPS) / release_interval);
    }
    else {
      // however, there must be at least one release
      return ceil_f(train_duration / release_interval);
    }
  }

  std::string release_pattern_to_data_model(Json::Value& mcell_node) const;
};

} // namespace mcell


#endif // SRC4_RELEASE_EVENT_H_
