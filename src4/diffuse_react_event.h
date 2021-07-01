/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_DIFFUSE_REACT_EVENT_H_
#define SRC4_DIFFUSE_REACT_EVENT_H_

#include <vector>
#include <deque>
#include "../libs/boost/container/small_vector.hpp"

#include "base_event.h"
#include "partition.h"
#include "collision_structs.h"


namespace MCell {

class World;
class Partition;
class Molecule;

// we are executing diffusion every iteration, do not change this constant
const double DIFFUSE_REACT_EVENT_PERIODICITY = 1.0;

// this is an arbitrary value but let's say that we do now want to
// simulate more than 100 iterations at once, this is used in scheduler
// to find the next barrier
const double DIFFUSION_TIME_UPPER_LIMIT = 100.0;

enum class RayTraceState {
  UNDEFINED,
  HIT_SUBPARTITION,
  RAY_TRACE_HIT_WALL,
  FINISHED
};

enum class WallRxnResult {
  INVALID,
  TRANSPARENT,
  REFLECT,
  DESTROYED
};


class TileNeighborVector: public std::deque<WallTileIndexPair> {
public:
  void dump(const std::string extra_comment, const std::string ind) {
    std::cout << ind << extra_comment << "\n";
    for (uint i = 0; i < size(); i++) {
      std::cout << ind << i << ": " << at(i).tile_index << "\n";
    }
  }
};


/**
 * Used as a pair molecule id, remaining timestep for molecules newly created in diffusion.
 * Using name action instead of event because events are handled by scheduler and are ordered by time.
 * These actions are simply processes in a queue (FIFO) manner.
 *
 * Used in diffuse_react _event_t and in partition_t.
 */
class DiffuseAction {
public:
  // position where this mol was created is
  // used to avoid rebinding for surf+vol->surf+vol reactions
  DiffuseAction(
      const molecule_id_t id_,
      const WallTileIndexPair& where_created_this_iteration_)
    :
      id(id_),
      where_created_this_iteration(where_created_this_iteration_) {
  }

  // position where the molecule was created may be unset when it was not a result of surface reaction
  DiffuseAction(const molecule_id_t id_)
    : id(id_) {
  }

  // defined because of usage in calendar_t
  const DiffuseAction& operator->() const { // TODO: remove
     return *this;
  }

  molecule_id_t id;

  // used to avoid rebinding for surf+vol->surf+vol reactions
  // TODO: can be removed and replaced with Molecule::previous_wall_index?
  WallTileIndexPair where_created_this_iteration;
};


/**
 * Diffuse all molecules with a given time step.
 * When a molecule is diffused, it is checked for collisions and reactions
 * are evaluated. Molecules newly created in reactions and diffused with their
 * remaining time.
 */
class DiffuseReactEvent : public BaseEvent {
public:
  DiffuseReactEvent(World* world_) :
    BaseEvent(EVENT_TYPE_INDEX_DIFFUSE_REACT),
    world(world_), time_up_to_next_barrier(FLT_INVALID) {

    // repeat this event each iteration
    periodicity_interval = DIFFUSE_REACT_EVENT_PERIODICITY;
  }

  void step() override;
  void dump(const std::string ind = "") const override;

  bool update_event_time_for_next_scheduled_time() override {
    assert(time_up_to_next_barrier != FLT_INVALID);
    // the next time to schedule
    if (time_up_to_next_barrier < periodicity_interval) {
      event_time = event_time + time_up_to_next_barrier;
    }
    else {
      event_time = event_time + periodicity_interval;
    }
    return true;
  }

  bool may_be_blocked_by_barrier_and_needs_set_time_step() const override {
    // DiffuseReactEvent must execute only up to a barrier such as CountEvent
    return true;
  }

  double get_max_time_up_to_next_barrier() const override {
    return DIFFUSION_TIME_UPPER_LIMIT;
  }

  void set_barrier_time_for_next_execution(const double time_up_to_next_barrier_) override {
    // scheduler says to this event for how long it can execute
    // either the maximum time step (periodicity_interval) or time up to the
    // first barrier
    release_assert(time_up_to_next_barrier_ > 0 && "Diffusion must advance even if a little bit");
    assert(cmp_eq(time_up_to_next_barrier_, round_f(time_up_to_next_barrier_)) &&
        "Time up to the next barrier is expected to be a whole number");
    time_up_to_next_barrier = time_up_to_next_barrier_;
  }

  void add_diffuse_action(const DiffuseAction& action) {
    new_diffuse_actions.push_back(action);
  }

  bool before_this_iterations_end(const double time) const {
    return cmp_lt(time, event_time + periodicity_interval, EPS);
  }

  // used also directly from MCell API
  // returns true if molecule survived
  bool outcome_unimolecular(
      Partition& p,
      Molecule& vm,
      const double scheduled_time,
      BNG::RxnClass* rxn_class,
      const BNG::rxn_class_pathway_index_t pathway_index,
      MoleculeIdsVector* optional_product_ids = nullptr
  );

  // returns true if molecule survived
  bool cross_transparent_wall(
      Partition& p,
      const Collision& collision,
      Molecule& vm, // moves vm to the reflection point
      Vec3& remaining_displacement,
      double& t_steps,
      double& elapsed_molecule_time,
      wall_index_t& last_hit_wall_index
  );

  World* world;

  // this event diffuses all molecules that have this diffusion time_step
  double time_up_to_next_barrier;

private:
  // auxiliary array used to store result from Partition::get_molecules_ready_for_diffusion
  // using the same array every iteration in order not to reallocate it every iteration
  MoleculeIdsVector molecules_ready_array;

  // internal event's schedule of molecules newly created in reactions that must be diffused
  std::vector<DiffuseAction> new_diffuse_actions;

  double get_max_time(Partition& p, Molecule& m);

  void diffuse_molecules(Partition& p, const MoleculeIdsVector& indices);

  void diffuse_single_molecule(
      Partition& p,
      const molecule_id_t vm_id,
      WallTileIndexPair where_created_this_iteration
  );

  // ---------------------------------- volume molecules ----------------------------------
  void diffuse_vol_molecule(
      Partition& p,
      Molecule& vm,
      double& max_time,
      const double diffusion_start_time,
      WallTileIndexPair& where_created_this_iteration
  );

  bool collide_and_react_with_vol_mol(
      Partition& p,
      const Collision& collision,
      Vec3& displacement,
      const double t_steps,
      const double r_rate_factor,
      const double elapsed_molecule_time
  );

  int collide_and_react_with_surf_mol(
      Partition& p,
      const Collision& collision,
      const double r_rate_factor,
      WallTileIndexPair& where_created_this_iteration,
      wall_index_t& last_hit_wall_index,
      Vec3& remaining_displacement,
      double& t_steps,
      double& elapsed_molecule_time
  );

  WallRxnResult collide_and_react_with_walls(
      Partition& p,
      Collision& collision,
      const double r_rate_factor,
      const double elapsed_molecule_time,
      const double t_steps
  );

  // ---------------------------------- surface molecules ----------------------------------

  void diffuse_surf_molecule(
      Partition& p,
      Molecule& sm,
      double& max_time,
      const double diffusion_start_time
  );

  wall_index_t ray_trace_surf(
      Partition& p,
      const BNG::Species& species,
      const molecule_id_t sm_id,
      Vec2& remaining_displacement,
      Vec2& new_pos,
      BNG::RxnClass*& absorb_now_rxn_class // set to non-null is molecule has to be absorbed
  );

  bool react_2D_all_neighbors(
      Partition& p,
      Molecule& sm,
      const double time, // same argument as t passed in mcell3 (come up with a better name)
      const double diffusion_start_time // diffusion_start_time + elapsed_molecule_time should be the time when reaction occurred
  );

  bool react_2D_intermembrane(
      Partition& p, Molecule& sm,
      const double t_steps, const double diffusion_start_time
  );

  // ---------------------------------- reactions ----------------------------------
  int find_surf_product_positions(
      Partition& p,
      const Collision& collision,
      const BNG::RxnRule* rxn,
      const Molecule* reacA, const bool keep_reacA,
      const Molecule* reacB, const bool keep_reacB,
      const Molecule* surf_reac,
      const BNG::RxnProductsVector& actual_products,
      GridPosVector& assigned_surf_product_positions,
      uint& num_surface_products,
      bool& surf_pos_reacA_is_used
  );

  int outcome_bimolecular(
      Partition& p,
      const Collision& collision,
      const int path,
      const double remaining_time_step
  );

  int outcome_intersect(
      Partition& p,
      BNG::RxnClass* rxn_class,
      const BNG::rxn_class_pathway_index_t pathway_index,
      Collision& collision,
      const double time
  );


  void handle_rxn_callback(
      Partition& p,
      const Collision& collision,
      const double time,
      const BNG::RxnRule* rxn,
      const Molecule* reac1,
      const Molecule* reac2,
      const MoleculeIdsVector& product_ids
  );

  orientation_t determine_orientation_depending_on_surf_comp(
      const species_id_t prod_species_id,
      const Molecule* surf_reac
  );

  int outcome_products_random(
      Partition& p,
      const Collision& collision,
      const double remaining_time_step,
      const BNG::rxn_class_pathway_index_t pathway_index,
      bool& keep_reacA,
      bool& keep_reacB,
      MoleculeIdsVector* optional_product_ids = nullptr
  );

  void pick_unimol_rxn_class_and_set_rxn_time(
      Partition& p,
      const double remaining_time_step,
      Molecule& vm
  );

  bool react_unimol_single_molecule(
      Partition& p,
      const molecule_id_t vm_id
  );
};


RayTraceState ray_trace_vol(
    Partition& p,
    rng_state& rng,
    const molecule_id_t vm_id, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const bool can_vol_react,
    const wall_index_t previous_reflected_wall, // is WALL_INDEX_INVALID when our molecule did not replect from anything this iddfusion step yet
    Vec3& remaining_displacement, // in/out - recomputed if there was a reflection
    CollisionsVector& molecule_collisions // possible reactions in this part of way marching, ordered by time
);

void sort_collisions_by_time(CollisionsVector& molecule_collisions);

} // namespace mcell

#endif // SRC4_DIFFUSE_REACT_EVENT_H_
