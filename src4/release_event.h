/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_RELEASE_EVENT_H_
#define SRC4_RELEASE_EVENT_H_

#include <vector>

#include "base_event.h"
#include "region_expr.h"

namespace Json {
class Value;
}

namespace MCell {

class Partition;
class Region;
class Wall;
class Grid;
class InitialSurfaceReleases;
class DiffuseReactEvent;

// TODO: replace with defines from API
enum class ReleaseShape {
  UNDEFINED = -1,  /* Not specified */
  SPHERICAL,       /* Volume enclosed by a sphere */
  // CUBIC,           /* Volume enclosed by a cube */ (might be supported, needs to be tested)
  // ELLIPTIC,        /* Volume enclosed by an ellipsoid */ (might be supported, needs to be tested)
  // RECTANGULAR,     /* Volume enclosed by a rect. solid */ (might be supported, needs to be tested)
  SPHERICAL_SHELL, /* Surface of a sphere */ // not tested yet
  REGION,          /* Inside/on the surface of an arbitrary region */

  // Individual mol. placement by list
  LIST,

  // Special case for MCell3 compatibility, supports handling of MOLECULE_NUMBER and MOLECULE_DENSITY
  // Same functionality as REGION, however implemented slightly differently
  INITIAL_SURF_REGION,
};


enum class ReleaseNumberMethod {
  INVALID,
  CONST_NUM, // used also for ReleaseShape::LIST
  GAUSS_NUM, // not supported yet
  VOL_NUM,
  CONCENTRATION_NUM,
  DENSITY_NUM
};

const int NUMBER_OF_TRAINS_UNLIMITED = -1;

// Data structure used to store LIST releases
class SingleMoleculeReleaseInfo {
public:
  SingleMoleculeReleaseInfo()
    : species_id(SPECIES_ID_INVALID), orientation(ORIENTATION_NONE),
      pos(POS_INVALID) {
  }

  species_id_t species_id;
  orientation_t orientation;
  Vec3 pos;
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
    release_number_method(ReleaseNumberMethod::INVALID),
    release_number(UINT_INVALID),
    concentration(FLT_INVALID),
    orientation(ORIENTATION_NONE),
    release_shape(ReleaseShape::UNDEFINED),
    location(FLT_INVALID),
    diameter(FLT_INVALID),
    region_llf(FLT_INVALID),
    region_urb(FLT_INVALID),
    // the default values for release pattern are such that there is
    // just one release
    delay(0),
    number_of_trains(1),
    train_interval(EPS),
    train_duration(EPS),
    release_interval(DBL_GIGANTIC),

    release_probability(1),

    actual_release_time(TIME_INVALID),
    current_train_from_0(0),
    current_release_in_train_from_0(0),
    world(world_),
    running_diffuse_event_to_update(nullptr) {
  }

  void step() override;

  // argument may be nullptr
  // NOTE: will need extra care for parallel diffusion
  void release_immediatelly(DiffuseReactEvent* running_diffuse_event_to_update_) {
      update_event_time_for_next_scheduled_time();
      running_diffuse_event_to_update = running_diffuse_event_to_update_;
      step();
  }

  // release events must be sorted by the actual release time as well
  bool needs_secondary_ordering() override {
    return true;
  }

  double get_secondary_ordering_value() override {
    assert(actual_release_time != TIME_INVALID);
    return actual_release_time;
  }

  // handles release patterns,
  // WARNING: must be called before the first insertion into the scheduler,
  // when no release pattern is set, simply sets event_time to 0
  // and on the second call it returns false
  bool update_event_time_for_next_scheduled_time() override;

  // DiffuseReactEvent must execute only up to this event
  // for MCell3 compatibility but otherwise not sure why this is necessary
  bool is_barrier() const override { return true; }

  void dump(const std::string indent) const override;
  void to_data_model(Json::Value& mcell_node) const override;

  bool needs_release_pattern() const {
    bool single_release_at_t0 = delay == 0 && number_of_trains == 1 && get_num_releases_per_train() == 1;
    return !single_release_at_t0;
  }

public:
  std::string release_site_name; // name of releaser site from which was this event created

  species_id_t species_id;

  ReleaseNumberMethod release_number_method; // specifies what does the release_number mean
  uint release_number; // number of molecules to release
  double concentration; // or density for surface releases

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
  // defines walls of a region for surface release
  // TODO: replace with some a better structure
  std::vector<CummAreaPWallIndexPair> cumm_area_and_pwall_index_pairs;

  RegionExpr region_expr;

  // for volume molecule releases into a region
  Vec3 region_llf; // note: this is fully specified by the region above, maybe remove in the future
  Vec3 region_urb; // note: this is fully specified by the region above as well

  // used when release_shape is ReleaseShape::List
  std::vector<SingleMoleculeReleaseInfo> molecule_list;

  std::string release_pattern_name;

  // --- release pattern information ---
  double delay;
  int number_of_trains; // -1 means that the number of trains is unlimited
  double train_interval;
  double train_duration;
  double release_interval;

  // This release does not occur every time, but rather with probability p.
  // Either the whole release occurs or none of it does; the probability does not
  // apply molecule-by-molecule. release_probability must be in the interval [0, 1].
  double release_probability;

private:
  double actual_release_time;

  // both values are initialized to 0 and counted from 0
  int current_train_from_0;
  int current_release_in_train_from_0;

  World* world;

  // if not nullptr, we need to inform this event that there are new
  // molecules to be released
  DiffuseReactEvent* running_diffuse_event_to_update;

private:
  bool skip_due_to_release_probability();

  // may end with error
  void report_release_failure(const std::string& msg);

  uint calculate_number_to_release();

  int randomly_remove_molecules(
      Partition& p, const MoleculeIdsVector& mol_ids_in_region, int number_to_remove);

  // for surface molecule releases
  int vacuum_from_regions(int number_to_remove);
  void release_onto_regions(int& computed_release_number);

  // for volume molecule releases into a region
  bool is_point_inside_region_expr_recursively(
      Partition& p, const Vec3& pos, const RegionExprNode* region_expr_node
  );
  uint num_vol_mols_from_conc(bool &exact_number);
  int vacuum_inside_regions(int number_to_remove);
  void release_inside_regions(int& computed_release_number);

  // for volume molecule releases
  void release_ellipsoid_or_rectcuboid(int computed_release_number);

  // for list releases
  void release_list();

  // for releases specified by MODIFY_SURFACE_REGIONS -> MOLECULE_NUMBER or MOLECULE_DENSITY
  void init_surf_mols_by_number(
      Partition& p, const Region& reg, const InitialSurfaceReleases& info);
  void init_surf_mols_by_density(
      Partition& p, Wall& w, std::map<species_id_t, uint>& num_released_per_species);
  void release_initial_molecules_onto_surf_regions();

  void schedule_for_immediate_diffusion_if_needed(
      const molecule_id_t id, const WallTileIndexPair& where_released = WallTileIndexPair());

  double get_release_delay_time() const {
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

  void to_data_model_as_one_release_site(
      Json::Value& mcell_node,
      const species_id_t species_id_override,
      const orientation_t orientation_override,
      // points_list indices are set only when
      // release_shape == ReleaseShape::LIST
      const std::string& name_override,
      const uint points_list_begin_index,
      const uint points_list_end_index
  ) const;
};


// utilities used also from ConcentrationClampReleaseEvent
size_t cum_area_bisect_high(const std::vector<CummAreaPWallIndexPair>& array, double val);
void dump_cumm_area_and_pwall_index_pairs(
    const std::vector<CummAreaPWallIndexPair>& cumm_area_and_pwall_index_pairs, const std::string ind);


} // namespace mcell


#endif // SRC4_RELEASE_EVENT_H_
