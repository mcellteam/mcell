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

#ifndef SRC4_RAY_TRACER_H_
#define SRC4_RAY_TRACER_H_

#include "/home/ahusar/src4_ref/embree_test/embree_debug_install/include/embree3/rtcore.h"

#include "collision_structs.h"

struct rng_state;

namespace MCell {

class Partition;
class Molecule;

struct RayTraceMoleculeData {
  molecule_id_t molecule_id;
  species_id_t species_id;
};


struct HitInfo {
  uint molecule_id; // molecule ID that was hit
  float ray_tfar;
  // probably additional info is needed
};


class RayTracer {
public:
  RayTracer(Partition& p_)
    : p(p_), initialized(false), device(nullptr), scene(nullptr), walls_created(false), walls_geometry_id(0) {
  }
  ~RayTracer();

  // for now the geometry is fixed
  void initialize_and_create_geometry();

  void add_molecule(Molecule& vm);
  void update_molecule_position(Molecule& vm);
  void remove_molecule(Molecule& vm);

  RayTraceState ray_trace_vol(
      rng_state& rng,
      const molecule_id_t vm_id, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
      const wall_index_t last_hit_wall_index, // is WALL_INDEX_INVALID when our molecule did not reflect from anything this diffusion step yet
      Vec3& remaining_displacement, // in/out - recomputed if there was a reflection
      collision_vector_t& collisions // both mol mol and wall collisions
  );

private:
  void store_molecule_collision(
      const HitInfo& h, const Molecule& diffused_vm, const Vec3& displacement, collision_vector_t& collisions
  );

  void store_wall_collision(
      const wall_index_t wall_index, const Molecule& diffused_vm, Vec3& displacement,
      collision_vector_t&collisions
  );

  Partition& p;

  bool initialized;
  RTCDevice device;
  RTCScene scene;

  bool walls_created;
  uint walls_geometry_id;

  // extra data for each molecule, managed by this class
  std::map<molecule_id_t, RayTraceMoleculeData> molecule_data;
};

} // namespace MCell

#endif /* SRC4_RAY_TRACER_H_ */
