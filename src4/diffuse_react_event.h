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

#ifndef SRC4_DIFFUSE_REACT_EVENT_H_
#define SRC4_DIFFUSE_REACT_EVENT_H_

#include <vector>
#include <deque>
#include "../libs/boost/container/small_vector.hpp"

#include "base_event.h"
#include "partition.h"
#include "collision_structs.h"

#define TEST 1

namespace MCell {

class World;
class Partition;
class Molecule;


enum class RayTraceState {
  UNDEFINED,
  HIT_SUBPARTITION,
  RAY_TRACE_HIT_WALL,
  FINISHED
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

typedef small_vector<wall_index_t> wall_indices_t;


/**
 * Used as a pair molecule id, remaining timestep for molecules newly created in diffusion.
 * Using name action instead of event because events are handled by scheduler and are ordered by time.
 * These actions are simply processes in a queue (FIFO) manner.
 *
 * Used in diffuse_react _event_t and in partition_t.
 */
class DiffuseOrUnimolRxnAction {
public:
  enum class Type {
    DIFFUSE,
    UNIMOL_REACT
  };

  // DIFFUSE action
  DiffuseOrUnimolRxnAction(
      const DiffuseOrUnimolRxnAction::Type type_,
      const molecule_id_t id_,
      const float_t scheduled_time_,
      const WallTileIndexPair& where_created_this_iteration_)
    :
      id(id_),
      scheduled_time(scheduled_time_),
      type(type_),
      unimol_rx(nullptr),
      where_created_this_iteration(where_created_this_iteration_) {

    assert(scheduled_time >= 0.0);
    assert(type == Type::DIFFUSE);
    // position where the molecule was created may be invalid when it was not a result of surface reaction
  }

  // UNIMOL_REACT action
  DiffuseOrUnimolRxnAction(
      const DiffuseOrUnimolRxnAction::Type type_,
      const molecule_id_t id_,
      const float_t scheduled_time_,
      const BNG::RxnClass* unimol_rx_)
    :
      id(id_),
      scheduled_time(scheduled_time_),
      type(type_),
      unimol_rx(unimol_rx_) {
    assert(scheduled_time >= 0.0);
    assert(type == Type::UNIMOL_REACT);
    assert(unimol_rx != nullptr);
  }

  // defined because of usage in calendar_t
  const DiffuseOrUnimolRxnAction& operator->() const {
     return *this;
  }

  molecule_id_t id;
  float_t scheduled_time; // this is the scheduled time
  Type type;

  // when type is UNIMOL_REACT
  const BNG::RxnClass* unimol_rx;

  // when type is DIFFUSE
  // used to avoid rebinding for surf+vol->surf+vol reactions
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
  DiffuseReactEvent(World* world_, const float_t diffusion_time_step_) :
    BaseEvent(EVENT_TYPE_INDEX_DIFFUSE_REACT),
    world(world_),
    diffusion_time_step(diffusion_time_step_) {

    // repeat this event each time step
    periodicity_interval = diffusion_time_step;
  }
  void step();
  void dump(const std::string indent);

  World* world;

  // this event diffuses all molecules that have this diffusion time_step
  float_t diffusion_time_step;

private:
  // molecules newly created in reactions
  std::vector<DiffuseOrUnimolRxnAction> new_diffuse_or_unimol_react_actions;

  void diffuse_molecules(Partition& p, const std::vector< molecule_id_t >& indices);

  void diffuse_single_molecule(
      Partition& p,
      const molecule_id_t vm_id,
      const float_t diffusion_start_time,
      WallTileIndexPair where_created_this_iteration
  );

  // ---------------------------------- volume molecules ----------------------------------
  // FIXME: unify argument names
  void diffuse_vol_molecule(
      Partition& p,
      const molecule_id_t vm_id,
      const float_t remaining_time_step,
      const float_t diffusion_start_time,
      WallTileIndexPair& where_created_this_iteration
  );

  bool collide_and_react_with_vol_mol(
      Partition& p,
      const Collision& collision,
      Vec3& displacement,
      const float_t remaining_time_step,
      const float_t r_rate_factor,
      const float_t molecule_time
  );

  int collide_and_react_with_surf_mol(
      Partition& p,
      const Collision& collision,
      const float_t remaining_time_step,
      const float_t r_rate_factor,
      const float_t current_molecule_time,
      WallTileIndexPair& where_created_this_iteration
  );

  // ---------------------------------- surface molecules ----------------------------------

  void diffuse_surf_molecule(
      Partition& p,
      const molecule_id_t sm_id,
      const float_t current_time,
      const float_t remaining_time_step
  );

  wall_index_t ray_trace_surf(
      Partition& p,
      const BNG::Species& species,
      Molecule& sm,
      Vec2& remaining_displacement,
      Vec2& new_pos
  );

  bool react_2D_all_neighbors(
      Partition& p,
      Molecule& sm,
      const float_t time, // same argument as t passed in mcell3 (come up with a better name)
      const float_t diffusion_start_time, // diffusion_start_time + elapsed_molecule_time should be the time when reaction occurred
      const float_t elapsed_molecule_time
  );

  // ---------------------------------- reactions ----------------------------------
  int find_surf_product_positions(
      Partition& p,
      const Molecule* reacA, const bool keep_reacA,
      const Molecule* reacB, const bool keep_reacB,
      const Molecule* surf_reac,
      const BNG::RxnRule* rxn,
      small_vector<GridPos>& assigned_surf_product_positions
  );

  int outcome_bimolecular(
      Partition& p,
      const Collision& collision,
      int path,
      float_t remaining_time_step
  );

	// returns true if molecule urvived
  bool outcome_unimolecular(
      Partition& p,
      Molecule& vm,
      const float_t scheduled_time,
      const BNG::RxnRule* unimol_rx
  );

  int outcome_products_random(
      Partition& p,
      const Collision& collision,
      const float_t remaining_time_step,
      const reaction_index_t reaction_index,
      bool& keep_reacA,
      bool& keep_reacB
  );

  void create_unimol_rx_action(
      Partition& p,
      Molecule& vm,
      float_t remaining_time_step
  );

  bool react_unimol_single_molecule(
      Partition& p,
      const molecule_id_t vm_id,
      const float_t scheduled_time,
      const BNG::RxnClass* unimol_rx
  );
};


//FIXME: move to some utility namespace
RayTraceState ray_trace_vol(
    Partition& p,
    rng_state& rng,
    Molecule& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const wall_index_t previous_reflected_wall, // is WALL_INDEX_INVALID when our molecule did not replect from anything this iddfusion step yet
    Vec3& remaining_displacement, // in/out - recomputed if there was a reflection
    collision_vector_t& molecule_collisions, // possible reactions in this part of way marching, ordered by time
    Vec3& new_position,
    subpart_index_t& new_subpartition_index
);

void sort_collisions_by_time(collision_vector_t& molecule_collisions);


} // namespace mcell

#endif // SRC4_DIFFUSE_REACT_EVENT_H_
