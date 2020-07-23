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

#include "xmmintrin.h"
#include "pmmintrin.h"


#include "ray_tracer.h"
#include "partition.h"

#include "logging.h"

#include "geometry_utils.h"
#include "geometry_utils.inc"
#include "collision_utils.inc"
#include "wall_utils.inc"

namespace MCell {

static void print_embree_error(const RTCError code)
{
  if (code == RTC_ERROR_NONE) {
    return;
  }

  cout << "Embree error code:";
  switch (code) {
    case RTC_ERROR_UNKNOWN          : cout << "RTC_ERROR_UNKNOWN"; break;
    case RTC_ERROR_INVALID_ARGUMENT : cout << "RTC_ERROR_INVALID_ARGUMENT"; break;
    case RTC_ERROR_INVALID_OPERATION: cout << "RTC_ERROR_INVALID_OPERATION"; break;
    case RTC_ERROR_OUT_OF_MEMORY    : cout << "RTC_ERROR_OUT_OF_MEMORY"; break;
    case RTC_ERROR_UNSUPPORTED_CPU  : cout << "RTC_ERROR_UNSUPPORTED_CPU"; break;
    case RTC_ERROR_CANCELLED        : cout << "RTC_ERROR_CANCELLED"; break;
    default                         : cout << "invalid error code"; break;
  }
  cout << "\n";
}


#define CHECK() do { \
    RTCError code = rtcGetDeviceError(device); \
    if (code != RTC_ERROR_NONE) { \
      print_embree_error(code); \
      assert(code == RTC_ERROR_NONE); \
      exit(1); \
    } \
  } while (0)


//#define CHECK() assert(rtcGetDeviceError(device) == RTC_ERROR_NONE);

// extra data used and created when intersecting hits
struct UserContext {
  uint diffused_species_id;

  uint last_hit_wall;

  void add_hit(const HitInfo& hit) {
    hits.push_back(hit);
  }

  std::vector<HitInfo> hits;
};


// intersection context passed to intersect calls
struct IntersectContext
{
  RTCIntersectContext rtc_context;
  mutable UserContext* user_context;
};


static void error_function(void* userPtr, enum RTCError error, const char* str) {
  mcell_log("embree error %d: %s", error, str);
}


static RTCDevice initialize_device() {

  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);


  RTCDevice device = rtcNewDevice("threads=1"); // disabling verbose causes an issue
  CHECK();
  if (!device) {
    mcell_log("error %d: cannot create device", rtcGetDeviceError(NULL));
  }

  rtcSetDeviceErrorFunction(device, error_function, NULL);
  CHECK();
  return device;
}


RTCScene initialize_scene(RTCDevice device)
{
  release_assert(sizeof(float) == 1*4);

  RTCScene scene = rtcNewScene(device);
  CHECK();
  rtcSetSceneFlags(scene,RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
  CHECK();
  return scene;
}


RayTracer::~RayTracer() {
  if (initialized) {
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);
  }
}


static void wall_intersection_filter(const RTCFilterFunctionNArguments* args) {
  assert(args->N == 1);

  UserContext* user_context = ((const IntersectContext*)args->context)->user_context;

  if (args->valid[0] != -1) {
    assert(false && "Invalid ray");
    return; // this is the default handling in Embree tutorials
  }

  RTCHit* hit = (RTCHit*)args->hit;
  if (hit->primID == user_context->last_hit_wall) {
    // ignore this wall because we hit it before
    args->valid[0] = 0;
  }
}


void RayTracer::initialize_and_create_geometry() {

  assert(!initialized);
  device = initialize_device();
  scene = initialize_scene(device);
  initialized = true;

  // represent each wall in the embree scene,
  // ordering of walls and indices is exactly the same as in the partition
  // we do not distinguish types of walls

  RTCGeometry walls_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
  CHECK();

  const std::vector<Vec3>& partition_verts = p.get_geometry_vertices();

  if (!partition_verts.empty()) {
    // it is not allowed to allocate a buffer of size 0

    // vertices
    float* vertices = (float*)rtcSetNewGeometryBuffer(
        walls_geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3*sizeof(float), // stride
        partition_verts.size() // count of vertices
    );
    CHECK();

    for (size_t i = 0; i < partition_verts.size(); i++) {
      const Vec3& v = partition_verts[i];
      size_t base = i*3;
      vertices[base + 0] = v.x;
      vertices[base + 1] = v.y;
      vertices[base + 2] = v.z;
    }

    // walls
    const std::vector<Wall>& partition_walls = p.get_walls();
    uint* indices = (uint*)rtcSetNewGeometryBuffer(
        walls_geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        3*sizeof(uint), // stride
        partition_walls.size() // count of walls
    );
    CHECK();

    for (size_t i = 0; i < partition_walls.size(); i++) {
      const Wall& w = partition_walls[i];
      size_t base = i*3;
      indices[base + 0] = w.vertex_indices[0];
      indices[base + 1] = w.vertex_indices[1];
      indices[base + 2] = w.vertex_indices[2];
    }

    // define a callback that skips hits of wall that we already hit
    rtcSetGeometryIntersectFilterFunction(walls_geom, wall_intersection_filter);
    CHECK();

    // we must commit geometry objects when they are set up
    rtcCommitGeometry(walls_geom);
    CHECK();

    // In rtcAttachGeometry(...), the scene takes ownership of the geom
    // by increasing its reference count. This means that we don't have
    // to hold on to the geom handle, and may release it. The geom object
    // will be released automatically when the scene is destroyed.
    //
    // rtcAttachGeometry() returns a geometry ID. We could use this to
    // identify intersected objects later on.
    walls_geometry_id = rtcAttachGeometry(scene, walls_geom);
    CHECK();
    rtcReleaseGeometry(walls_geom);
    CHECK();
    rtcCommitScene(scene);
    CHECK();
  }
  else {
    walls_geometry_id = UINT_INVALID;
  }
}



static void molecule_intersection_filter(const RTCFilterFunctionNArguments* args) {

  assert(args->N == 1);

  UserContext* user_context = ((const IntersectContext*)args->context)->user_context;

  if (args->valid[0] != -1) {
    assert(false && "Invalid ray");
    return; // this is the default handling in Embree tutorials
  }

  // we will ignore this hit
  args->valid[0] = 0;

  RTCRay* ray = (RTCRay*)args->ray;
  RTCHit* hit = (RTCHit*)args->hit;

  // aren't we over the set distance
  if (ray->tfar) {
    return;
  }

  // what species is this
  // TODO: , can we react?
  RayTraceMoleculeData* mols = (RayTraceMoleculeData*)args->geometryUserPtr;


  HitInfo hit_info;
  hit_info.molecule_id = mols[hit->primID].molecule_id;
  hit_info.ray_tfar = ray->tfar;
  user_context->add_hit(hit_info);

  /*
  printf("S - Found intersection on geometry %d, primitive %d at tfar=%f, species = %d, my_id = %d\n",
      hit->geomID,
      hit->primID,
      ray->tfar,
      user_context->diffused_species_id,
      mols[hit->primID].id
  );*/
}


void RayTracer::add_molecule(Molecule& vm) {
  assert(initialized);

  assert(vm.is_vol());

  RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);
  CHECK();

  float* sphere = (float*)rtcSetNewGeometryBuffer(
      geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4,
      4*sizeof(float), // stride
      1 // count
  );
  CHECK();

  sphere[0] = vm.v.pos.x;
  sphere[1] = vm.v.pos.y;
  sphere[2] = vm.v.pos.z;
  sphere[3] = p.config.rx_radius_3d; // the radius of the sphere is the reaction radius

  // set additional information for our sphere, only molecule id for now
  RayTraceMoleculeData data;
  data.molecule_id = vm.id;
  data.species_id = vm.species_id;
  // we need to maintain these data by ourselves,
  RayTraceMoleculeData& data_ref = (molecule_data[vm.id] = data);

  // make data available when calling intersection filter
  rtcSetGeometryUserData(geom, (void*)&data_ref);
  CHECK();

  // define a callback that stores information on molecule hits and continues with ray tracing
  rtcSetGeometryIntersectFilterFunction(geom, molecule_intersection_filter);
  CHECK();

  // and finally register our molecule to the scene
  rtcCommitGeometry(geom);
  CHECK();
  vm.v.embree_geom_id = rtcAttachGeometry(scene, geom);
  CHECK();
  rtcReleaseGeometry(geom);
  CHECK();
  rtcCommitScene(scene);
  CHECK();
}


void RayTracer::update_molecule_position(Molecule& vm) {
  assert(initialized);

  // get geometry assigned to this molecule
  assert(vm.v.embree_geom_id != UINT_INVALID);
  RTCGeometry geom = rtcGetGeometry(scene, vm.v.embree_geom_id);
  CHECK();
  assert(geom != nullptr);

  // update it (no need to change the radius)
  float* coords = (float*) rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0);
  CHECK();
  coords[0] = vm.v.pos.x;
  coords[1] = vm.v.pos.y;
  coords[2] = vm.v.pos.z;

  // and update the scene
  rtcUpdateGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0);
  CHECK();
  rtcCommitGeometry(geom);
  CHECK();
  rtcCommitScene(scene);
  CHECK();
}


void RayTracer::remove_molecule(Molecule& vm) {
  assert(initialized);

  rtcDetachGeometry(scene, vm.v.embree_geom_id);
  CHECK();
  rtcCommitScene(scene);
  CHECK();
}


void RayTracer::store_molecule_collision(
    const HitInfo& h, const Molecule& diffused_vm, const Vec3& displacement,
    collision_vector_t& collisions
) {

  const Molecule& colliding_vm = p.get_m(h.molecule_id);

  BNG::RxnClass* rxn_class =
      p.get_all_rxns().get_bimol_rxn_class(diffused_vm.species_id, colliding_vm.species_id);

  if (rxn_class == nullptr || rxn_class->get_num_reactions() == 0) {
    // no reaction possible
    return;
  }

  // recompute the hit using float_t
  // time in HitInfo is when the sphere was hit
  float_t time;
  Vec3 pos;
  bool hit_ok = CollisionUtil::collide_mol(diffused_vm, displacement, colliding_vm, p.config.rx_radius_3d, time, pos);
  assert(hit_ok);
  // TODO: some molecules are found even when collide_mol rejects them,
  // find out why
  if (!hit_ok) {
    return;
  }

  collisions.push_back(
      Collision(CollisionType::VOLMOL_VOLMOL, &p, diffused_vm.id, time, pos, colliding_vm.id, rxn_class)
  );
}


void RayTracer::store_wall_collision(
    const wall_index_t wall_index, const Molecule& diffused_vm, Vec3& displacement,
    collision_vector_t& collisions
) {

  const Wall& w = p.get_wall(wall_index);
  float_t time;
  Vec3 pos;

  CollisionType collision_type =
      CollisionUtil::collide_wall(p, diffused_vm.v.pos, w, p.aux_rng, true, displacement, time, pos);

  collisions.push_back(
      Collision(collision_type, &p, diffused_vm.id, time, pos, wall_index)
  );
}


RayTraceState RayTracer::ray_trace_vol(
    rng_state& rng,
    const molecule_id_t vm_id, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const wall_index_t last_hit_wall_index, // is WALL_INDEX_INVALID when our molecule did not reflect from anything this diffusion step yet
    Vec3& remaining_displacement, // in/out - recomputed if there was a reflection
    collision_vector_t& collisions // both mol mol and wall collisions
) {
  assert(initialized);
  collisions.clear();

  Molecule& vm = p.get_m(vm_id); // we need to update the position of this molecule

  IntersectContext context;
  UserContext user_context;
  // setup custom information
  user_context.diffused_species_id = vm.species_id;
  user_context.last_hit_wall = last_hit_wall_index;
  context.user_context = &user_context;
  rtcInitIntersectContext(&context.rtc_context);
  CHECK();

  RTCRayHit rayhit;
  rayhit.ray.org_x = vm.v.pos.x;
  rayhit.ray.org_y = vm.v.pos.y;
  rayhit.ray.org_z = vm.v.pos.z;
  rayhit.ray.dir_x = remaining_displacement.x;
  rayhit.ray.dir_y = remaining_displacement.y;
  rayhit.ray.dir_z = remaining_displacement.z;
  rayhit.ray.tnear = 0;
  // although we are setting here the max distance to be traveled,
  // we get back time [0..]
  float_t max_dist = len3(remaining_displacement) + EPS;
  rayhit.ray.tfar = max_dist;
  rayhit.ray.mask = -1;
  rayhit.ray.flags = 0;

  rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
  rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

  // TODO:
  // rtcSetGeometryMask

  // cast ray
  rtcIntersect1(scene, &context.rtc_context, &rayhit);
  CHECK();

  // handle molecule collisions
  uint_set<molecule_id_t> already_processed_mols;
  for (const HitInfo& h: user_context.hits) {
    assert(h.ray_tfar <= 1.0 && "Hits behind our displacement should have been already filtered out.");

    // we may also hit the front and the back side of a molecule sphere
    if (already_processed_mols.count(h.molecule_id) != 0) {
      continue;
    }
    store_molecule_collision(h, vm, remaining_displacement, collisions);
    already_processed_mols.insert(h.molecule_id);
  }

  // and check if we hit a wall
  RayTraceState res_state;
  if (rayhit.hit.geomID == walls_geometry_id &&
      rayhit.hit.primID != RTC_INVALID_GEOMETRY_ID &&
      rayhit.ray.tfar <= 1.0
  ) {
    store_wall_collision(rayhit.hit.primID, vm, remaining_displacement, collisions);

    res_state = RayTraceState::RAY_TRACE_HIT_WALL;
  }
  else {
    // no wall was hit
    vm.v.pos = vm.v.pos + remaining_displacement;
    vm.v.subpart_index = p.get_subpart_index(vm.v.pos);
    res_state = RayTraceState::FINISHED;
  }

  return res_state;
}


} // namespace MCell
