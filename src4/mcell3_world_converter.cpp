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

#include <stdarg.h>
#include <stdlib.h>
#include <set>

#include "bng/bng.h"

#include "logging.h"
#include "mcell_structs.h"
#include "mcell3_world_converter.h"

#include "world.h"
#include "release_event.h"
#include "diffuse_react_event.h"
#include "viz_output_event.h"
#include "mol_or_rxn_count_event.h"
#include "count_buffer.h"

#include "datamodel_defines.h"

#include "dump_state.h"

using namespace std;
using namespace BNG;

// checking major conversion blocks
#define CHECK(cond) do { if(!(cond)) { mcell_log_conv_error("Returning from %s after conversion error.\n", __FUNCTION__); return false; } } while (0)

// checking assumptions
#define CHECK_PROPERTY(cond) do { if(!(cond)) { mcell_log_conv_error("Expected '%s' is false. (%s - %s:%d)\n", #cond, __FUNCTION__, __FILE__, __LINE__); assert(false); return false; } } while (0)

// asserts - things that can never occur and will 'never' be supported


// holds global class
MCell::MCell3WorldConverter g_converter;


bool mcell4_convert_mcell3_volume(volume* s) {
  return g_converter.convert(s);
}


bool mcell4_run_simulation(const bool dump_initial_state) {
  g_converter.world->run_simulation(dump_initial_state);
  return true;
}

void mcell4_convert_to_datamodel() {
  g_converter.world->export_visualization_datamodel(DEFAULT_DATAMODEL_FILENAME);
  mcell_log("Datamodel was exported to '%s'.", DEFAULT_DATAMODEL_FILENAME);
}

void mcell4_delete_world() {
  g_converter.reset();
}


void mcell_log_conv_warning(char const *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  string fmt_w_warning = string ("Conversion warning: ") + fmt;
  mcell_logv_raw(fmt_w_warning.c_str(), args);
  va_end(args);
}


void mcell_log_conv_error(char const *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  string fmt_w_warning = string ("Conversion error: ") + fmt;
  mcell_logv_raw(fmt_w_warning.c_str(), args);
  va_end(args);
}


namespace MCell {

static const char* get_sym_name(const sym_entry *s) {
  assert(s != nullptr);
  assert(s->name != nullptr);
  return s->name;
}


static mat4x4 t_matrix_to_mat4x4(const double src[4][4]) {
  mat4x4 res;

  for (int x = 0; x < 4; x++) {
    for (int y = 0; y < 4; y++) {
      res[x][y] = src[x][y];
    }
  }

  return res;
}


void MCell3WorldConverter::reset() {
  delete world;
  world = nullptr;
  mcell3_species_id_map.clear();
}


bool MCell3WorldConverter::convert(volume* s) {

  world = new World();

  CHECK(convert_simulation_setup(s));

  CHECK(convert_species(s));
  create_diffusion_events(); // cannot fail?
  CHECK(convert_rxns(s));

  // at this point, we need to create the first (and for now the only) partition
  // create initial partition with center at 0,0,0 - we woud like to have the partitions all the same,
  // not depend on some random initialization
  partition_id_t index = world->add_partition(Vec3(0, 0, 0));
  assert(index == PARTITION_ID_INITIAL);

  // convert geometry already puts geometry objects into partitions
  CHECK(convert_geometry_objects(s));

  // release events require wall information
  CHECK(convert_release_events(s));
  CHECK(convert_viz_output_events(s));
  CHECK(convert_mol_count_and_rxn_count_events(s));

  return true;
}


static float_t get_largest_abs_value(const vector3& v) {
  float_t max = 0;
  if (fabs_f(v.x) > max) {
    max = fabs_f(v.x);
  }
  if (fabs_f(v.y) > max) {
    max = fabs_f(v.y);
  }
  if (fabs_f(v.z) > max) {
    max = fabs_f(v.z);
  }
  return max;
}


static float_t get_largest_distance_from_center(const vector3& llf, const vector3& urb) {
  float_t max1 = get_largest_abs_value(llf);
  float_t max2 = get_largest_abs_value(urb);
  return max1 > max2 ? max1 : max2;
}


static uint get_even_higher_or_same_value(const uint val) {
  if (val % 2 == 0) {
    return val;
  }
  else {
    return val + 1;
  }
}


static float_t get_partition_edge_length(const World* world, const float_t largest_mcell3_distance_from_center) {
  // some MCell models have their partition boundary set exactly. we need to add a bit of margin
  return (largest_mcell3_distance_from_center + PARTITION_EDGE_EXTRA_MARGIN_UM / world->config.length_unit) * 2 ;
}


/**
 * Partition conversion greatly simplifies the variability in MCell3 where the partition can be an arbitrary box.
 * Here, it must be a cube and the first partition must be placed the way so that the coordinate origin
 * is in its corner.
 * The assumption (quite possibly premature) is that there will be big systems that are simulated
 * and the whole space will be split into multiple partitions anyway. And we do not care about the
 * number of subpartitions, they should only take up memory, not computation time.
 */
bool MCell3WorldConverter::convert_simulation_setup(volume* s) {
  // TODO_CONVERSION: many items are not checked
  world->total_iterations = s->iterations;
  world->config.time_unit = s->time_unit;
  world->config.length_unit = s->length_unit;
  world->config.grid_density = s->grid_density;
  world->config.rx_radius_3d = s->rx_radius_3d;
  world->config.vacancy_search_dist2 = s->vacancy_search_dist2 / s->length_unit;
  world->seed_seq = s->seed_seq;
  world->rng = *s->rng;

  // there seems to be just one partition in MCell but we interpret it as mcell4 partition size
  if (s->partitions_initialized) {
    // assuming that the mcell's bounding box if it is bigger than the partition
    if (
        s->partition_llf.x <= s->bb_llf.x || s->bb_urb.x <= s->partition_urb.x ||
        s->partition_llf.y <= s->bb_llf.y || s->bb_urb.y <= s->partition_urb.y ||
        s->partition_llf.z <= s->bb_llf.z || s->bb_urb.z <= s->partition_urb.z
    ) {
      mcell_log("Partitioning was specified, but it is smaller than the automatically determined bounding box.");

      float_t lu = s->length_unit;
      mcell_log("Bounding box in microns: [ %f, %f, %f ], [ %f, %f, %f ]",
          s->bb_llf.x*lu, s->bb_llf.y*lu, s->bb_llf.z*lu, s->bb_urb.x*lu, s->bb_urb.y*lu, s->bb_urb.z*lu);
      mcell_log("Partition box in microns: [ %f, %f, %f ], [ %f, %f, %f ]",
          s->partition_llf.x*lu, s->partition_llf.y*lu, s->partition_llf.z*lu, s->partition_urb.x*lu, s->partition_urb.y*lu, s->partition_urb.z*lu);
    }

    CHECK_PROPERTY(s->partition_urb.x > s->partition_llf.x);
    CHECK_PROPERTY(s->partition_urb.y > s->partition_llf.y);
    CHECK_PROPERTY(s->partition_urb.z > s->partition_llf.z);

    // simply choose the largest abs value from all the values
    float_t largest_mcell3_distance_from_center = get_largest_distance_from_center(s->partition_llf, s->partition_urb);

    world->config.partition_edge_length = get_partition_edge_length(world, largest_mcell3_distance_from_center);

    mcell_log("Using manually specified partition size (with margin): %f.",
      (double)world->config.partition_edge_length * world->config.length_unit);

    // number of subparts must be even so that the central subparts are aligned with the axes and not shifted
    world->config.subpartitions_per_partition_dimension = get_even_higher_or_same_value(s->num_subparts);
  }
  else {
    CHECK_PROPERTY(s->bb_urb.x >= s->bb_llf.x);
    CHECK_PROPERTY(s->bb_urb.y >= s->bb_llf.y);
    CHECK_PROPERTY(s->bb_urb.z >= s->bb_llf.z);

    // use MCell's bounding box, however, we must make a cube out of it
    float_t largest_mcell3_distance_from_center = get_largest_distance_from_center(s->bb_llf, s->bb_urb);
    float_t auto_length = get_partition_edge_length(world, largest_mcell3_distance_from_center);

    if (auto_length > PARTITION_EDGE_LENGTH_DEFAULT_UM / world->config.length_unit) {
      world->config.partition_edge_length = auto_length;
      mcell_log("Automatically determined partition size: %f.",
        (double)auto_length * world->config.length_unit);
    }
    else {
      world->config.partition_edge_length = PARTITION_EDGE_LENGTH_DEFAULT_UM / world->config.length_unit;
      mcell_log("Automatically determined partition size %f is smaller than default %f, using default.",
        (double)auto_length * world->config.length_unit, PARTITION_EDGE_LENGTH_DEFAULT_UM);
    }

    // nx_parts counts the number of boundaries, not subvolumes, also, there are always 2 extra subvolumes on the sides in mcell3
    int max_n_p_parts = max3_i(IVec3(s->nx_parts, s->ny_parts, s->nz_parts));
    world->config.subpartitions_per_partition_dimension = get_even_higher_or_same_value(max_n_p_parts - 3);
  }

  float_t l = world->config.partition_edge_length / 2 * s->length_unit;
  mcell_log("MCell4 partition bounding box in microns: [ %f, %f, %f ], [ %f, %f, %f ]", -l, -l, -l, l, l, l);


  world->config.randomize_smol_pos = s->randomize_smol_pos;

  CHECK_PROPERTY(s->dynamic_geometry_molecule_placement == 0
      && "DYNAMIC_GEOMETRY_MOLECULE_PLACEMENT '=' NEAREST_TRIANGLE is not supported yet"
  );

  world->config.use_expanded_list = s->use_expanded_list;

	// compute other constants
  world->config.init();

  return true;
}


// "" as expected_name is that the name is not checked
bool check_meta_object(geom_object* o, string expected_name = "") {
  assert(o != nullptr);

  if (expected_name != "") {
    CHECK_PROPERTY(get_sym_name(o->sym) == expected_name);
  }
  else {
    // just try that the name is set in debug mode
    get_sym_name(o->sym);
  }
  // root->last_name - not checked, contains some nonsense anyway
  CHECK_PROPERTY(o->object_type == META_OBJ);
  CHECK_PROPERTY(o->contents == nullptr);
  CHECK_PROPERTY(o->num_regions == 0);
  CHECK_PROPERTY(o->regions == nullptr);
  CHECK_PROPERTY(o->walls == nullptr);
  CHECK_PROPERTY(o->wall_p == nullptr);
  CHECK_PROPERTY(o->vertices == nullptr);
  CHECK_PROPERTY(o->total_area == 0);
  CHECK_PROPERTY(o->n_tiles == 0);
  CHECK_PROPERTY(o->n_occupied_tiles == 0);
  CHECK_PROPERTY(o->n_occupied_tiles == 0);
  CHECK_PROPERTY(t_matrix_to_mat4x4(o->t_matrix) == mat4x4(1) && "only identity matrix for now");
  // root->is_closed - not checked
  CHECK_PROPERTY(o->periodic_x == false);
  CHECK_PROPERTY(o->periodic_y == false);
  CHECK_PROPERTY(o->periodic_z == false);
  return true;
}


bool MCell3WorldConverter::convert_geometry_objects(volume* s) {

  geom_object* root_instance = s->root_instance;
  CHECK_PROPERTY(check_meta_object(root_instance, "WORLD_INSTANCE"));
  CHECK_PROPERTY(root_instance->next == nullptr);

  for (geom_object* instantiate_obj = root_instance->first_child; instantiate_obj != nullptr; instantiate_obj = instantiate_obj->next) {
    CHECK_PROPERTY(check_meta_object(instantiate_obj));

    // this is the name of the INSTANTIATE section
    std::string instantiate_name = get_sym_name(instantiate_obj->sym);

    // walls reference each other, therefore we must first create
    // empty wall objects in partitions,
    geom_object* curr_obj = instantiate_obj->first_child;
    while (curr_obj != nullptr) {
      if (curr_obj->object_type == POLY_OBJ) {
        create_uninitialized_walls_for_polygonal_object(curr_obj);
      }

      curr_obj = curr_obj->next;
    }

    // once all wall were created and mapping established,
    // we can fill-in all objects
    curr_obj = instantiate_obj->first_child;
    while (curr_obj != nullptr) {
      if (curr_obj->object_type == POLY_OBJ) {
        CHECK(convert_polygonal_object(curr_obj, instantiate_name));
      }
      else if (curr_obj->object_type == REL_SITE_OBJ) {
        // ignored
      }
      else {
        CHECK_PROPERTY(false && "Unexpected type of object");
      }

      curr_obj = curr_obj->next;
    }

  } // for each scene/INSTANTIATE section

  return true;
}


// we do not check anything that might not be supported from the mcell3 side,
// the actual checks are in convert_polygonal_object
void MCell3WorldConverter::create_uninitialized_walls_for_polygonal_object(const geom_object* o) {
  GeometryObject obj;

  // create objects for each wall
  for (int i = 0; i < o->n_walls; i++) {
    wall* w = o->wall_p[i];

    // which partition?
    partition_id_t partition_id = world->get_partition_index(*w->vert[0]);
    if (partition_id == PARTITION_ID_INVALID) {
      Vec3 v(*w->vert[0]);
      v = v * Vec3(world->config.length_unit);
      mcell_error("Error: vertex %s does not fit any partition.", v.to_string().c_str());
    }

    // check that the remaining vertices are in the same partition
    for (uint k = 1; k < VERTICES_IN_TRIANGLE; k++) {
      partition_id_t curr_partition_index = world->get_partition_index(*w->vert[k]);

      if (partition_id != curr_partition_index) {
        Vec3 pos(*w->vert[k]);
        pos = pos * Vec3(world->config.length_unit);
        mcell_error("Error: whole walls must be in a single partition is for now, vertex %s is out of bounds", pos.to_string().c_str());
      }
    }

    // create the wall in that partition but do not set anything else yet
    Partition& p = world->get_partition(partition_id);
    Wall& new_wall = p.add_uninitialized_wall(world->get_next_wall_id());

    // remember mapping
    add_mcell4_wall_index_mapping(w, PartitionWallIndexPair(partition_id, new_wall.index));
  }
}


bool MCell3WorldConverter::convert_wall_and_update_regions(
    const wall* w, GeometryObject& object,
    const region_list* rl
) {

  PartitionWallIndexPair wall_pindex = get_mcell4_wall_index(w);
  Partition& p = world->get_partition(wall_pindex.first);
  Wall& wall = p.get_wall(wall_pindex.second);

  // bidirectional mapping
  wall.object_id = object.id;
  wall.object_index = object.index;
  object.wall_indices.push_back(wall.index);

  wall.side = w->side;

  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    // this vertex was inserted into the same partition as the whole object
    PartitionVertexIndexPair vert_pindex = get_mcell4_vertex_index(w->vert[i]);
    assert(wall_pindex.first == vert_pindex.first);
    wall.vertex_indices[i] = vert_pindex.second;
  }

  wall.uv_vert1_u = w->uv_vert1_u;
  wall.uv_vert2 = w->uv_vert2;
  wall.area = w->area;
  wall.normal = w->normal;
  wall.unit_u = w->unit_u;
  wall.unit_v = w->unit_v;
  wall.distance_to_origin = w->d;
  wall.wall_constants_precomputed = true;

  for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
    edge* e = w->edges[i];
    assert(e != nullptr);

    Edge& edge = wall.edges[i];

    if (e->forward != nullptr) {
      edge.forward_index = get_mcell4_wall_index(e->forward).second;
    }
    else {
      edge.forward_index = WALL_INDEX_INVALID;
    }
    if (e->backward != nullptr) {
      edge.backward_index = get_mcell4_wall_index(e->backward).second;
    }
    else {
      edge.backward_index = WALL_INDEX_INVALID;
    }
    edge.set_translate(e->translate);
    edge.set_cos_theta(e->cos_theta);
    edge.set_sin_theta(e->sin_theta);

    // e->edge_num_used_for_init is valid only if both fwd and backw walls are present
    if (e->forward != nullptr && e->backward != nullptr) {
      edge.edge_num_used_for_init = e->edge_num_used_for_init; // added only for mcell4
    }
    else {
      edge.edge_num_used_for_init = EDGE_INDEX_INVALID;
    }
  }

  for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
    if (w->nb_walls[i] != nullptr) {
#ifndef NDEBUG
      PartitionWallIndexPair pindex = get_mcell4_wall_index(w->nb_walls[i]);
      assert(pindex.first == wall_pindex.first && "Neighbors must be in the same partition for now");
#endif
      wall.nb_walls[i] = get_mcell4_wall_index(w->nb_walls[i]).second;
    }
    else {
      wall.nb_walls[i] = WALL_INDEX_INVALID;
    }
  }

  CHECK_PROPERTY(w->grid == nullptr); // for now

  // walls use some of the flags used by species
  CHECK_PROPERTY(
      w->flags == COUNT_CONTENTS ||
      w->flags == COUNT_RXNS ||
      w->flags == (COUNT_CONTENTS | COUNT_RXNS) ||
      w->flags == (COUNT_CONTENTS | COUNT_ENCLOSED) ||
      w->flags == (COUNT_RXNS | COUNT_ENCLOSED) ||
      w->flags == (COUNT_CONTENTS | COUNT_RXNS | COUNT_ENCLOSED)
  );

  // now let's handle regions
  std::set<species_id_t> region_species_from_mcell3;
  for (const region_list *r = rl; r != NULL; r = r->next) {
    region* reg = r->reg;
    assert(reg != nullptr);

    // are we in this region?
    if (get_bit(reg->membership, wall.side)) {
      auto pindex = get_mcell4_region_index(reg);
      CHECK_PROPERTY(pindex.first == wall_pindex.first);

      region_index_t region_index = pindex.second;
      wall.regions.insert_unique(region_index);

      // add our wall to the region
      Region& mcell4_reg = p.get_region(region_index);
      if (mcell4_reg.species_id != SPECIES_ID_INVALID) {
        region_species_from_mcell3.insert(mcell4_reg.species_id);
      }

      // which of our walls are region edges?
      set<edge_index_t> edge_indices;
      if (reg->boundaries != nullptr) {
        for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
          edge* e = w->edges[i];
          uint keyhash = (unsigned int)(intptr_t)(e);
          void* key = (void *)(e);
          if (pointer_hash_lookup(reg->boundaries, key, keyhash)) {
            edge_indices.insert(i);
          }
        }
      }

      assert(mcell4_reg.walls_and_edges.count(wall.index) == 0);
      mcell4_reg.walls_and_edges[wall.index] = edge_indices;
    }
  }

  // check that associated regions used the same 'surf_classes'/species
  std::set<species_id_t> wall_species_from_mcell3;
  for (surf_class_list *sl = w->surf_class_head; sl != nullptr; sl = sl->next) {
    CHECK_PROPERTY(sl->surf_class != nullptr);
    wall_species_from_mcell3.insert(get_mcell4_species_id(sl->surf_class->species_id));
  }
  CHECK_PROPERTY(wall_species_from_mcell3 == region_species_from_mcell3);

  // finally, we must let the partition know that
  // we initialized the wall
  p.finalize_wall_creation(wall.index);

  return true;
}

// result is not used anywhere
bool MCell3WorldConverter::convert_region(Partition& p, const region* r, region_index_t& region_index) {

  assert(r != nullptr);

  Region new_region;

  new_region.name = get_sym_name(r->sym);

  // u_int hashval;          // ignored
  // char *region_last_name; // ignored
  // struct geom_object *parent; // ignored
  CHECK_PROPERTY(r->element_list_head == nullptr); // should be used only at parse time

  // r->membership - Each bit indicates whether the corresponding wall is in the region */
  // and r->boundaries - edges of region are handled in convert_wall_and_update_regions

  CHECK_PROPERTY(r->sm_dat_head == nullptr); // should be null during initial conversion

  if (r->surf_class != nullptr) {
    new_region.species_id = get_mcell4_species_id(r->surf_class->species_id);
  }
  else {
    new_region.species_id = SPECIES_ID_INVALID;
  }

  CHECK_PROPERTY(r->region_has_all_elements || r->bbox == nullptr); // so far it can be set only to all-encopasing region
  // CHECK_PROPERTY(r->area == 0);  // ignored for now
  // CHECK_PROPERTY(r->volume == 0);  // ignored for now
  // CHECK_PROPERTY(r->flags == 0); // ignored for now
  // CHECK_PROPERTY(r->manifold_flag == 0); // ignored, do we care?

  string parent_name = get_sym_name(r->parent->sym);
  const GeometryObject* obj = world->get_partition(0).find_geometry_object_by_name(parent_name);
  CHECK_PROPERTY(obj != nullptr);
  new_region.geometry_object_id = obj->id;

  new_region.id = world->get_next_region_id();
  region_index = p.add_region_and_set_its_index(new_region);
  return true;
}


bool MCell3WorldConverter::convert_polygonal_object(const geom_object* o, const std::string& instantiate_name) {

  // --- object ---

  // we already checked in create_uninitialized_walls_for_polygonal_object
  // that the specific walls of this fit into a single partition
  // TODO_CONVERSION: improve this check for the whole object
  partition_id_t partition_index = world->get_partition_index(*o->vertices[0]);
  Partition& p = world->get_partition(partition_index);

  GeometryObject& obj = p.add_uninitialized_geometry_object(world->get_next_geometry_object_id());

  // o->next - ignored
  // o->parent - ignored
  CHECK_PROPERTY(o->first_child == nullptr);
  CHECK_PROPERTY(o->last_child == nullptr);
  obj.name = get_sym_name(o->sym);
  obj.parent_name = instantiate_name;
  // o->last_name - ignored
  CHECK_PROPERTY(o->object_type == POLY_OBJ);
  CHECK_PROPERTY(o->contents != nullptr); // ignored for now, not sure what is contained

  CHECK_PROPERTY(o->n_walls == o->n_walls_actual); // ignored
  CHECK_PROPERTY(o->walls == nullptr); // this is null for some reason
  CHECK_PROPERTY(o->wall_p != nullptr);

  // --- regions ---
  uint reg_cnt = 0;
  for (region_list *r = o->regions; r != NULL; r = r->next) {
    region* reg = r->reg;
    assert(reg != nullptr);

    region_index_t region_index;
    CHECK(convert_region(p, reg, region_index));

    add_mcell4_region_index_mapping(reg, PartitionRegionIndexPair(partition_index, region_index));
    reg_cnt++;
  }
  CHECK_PROPERTY(o->num_regions == reg_cnt);
  // TODO: add regions to the object

  // --- vertices ---
  // to stay identical to mcell3, will use the exact number of vertices as in mcell3, for this to work,
  // vector_ptr_to_vertex_index_map is a 'global' map for the whole conversion process
  // one of the reasons to not to copy vertex coordinates is that they are shared among triangles of an object
  // and when we move one vertex of the object, we transform all the triangles (walls) that use it
  for (int i = 0; i < o->n_verts; i++) {
    // insert vertex into the right partition and returns partition index and vertex index
    vertex_index_t new_vertex_index = p.add_geometry_vertex(*o->vertices[i]);
    add_mcell4_vertex_index_mapping(o->vertices[i], PartitionVertexIndexPair(partition_index, new_vertex_index));
  }

  // --- walls ---

  // vertex info contains also partition indices when it is inserted into the
  // world geometry

  //walls_vertices.resize(o->n_walls);
  for (int i = 0; i < o->n_walls; i++) {
    //walls_vertices[i].resize(VERTICES_IN_TRIANGLE);
    // uses precomputed map vector_ptr_to_vertex_index_map to transform vertices
    // also uses mcell3_region_to_mcell4_index mapping to set that it belongs to a given region
    //
    CHECK(convert_wall_and_update_regions(o->wall_p[i], obj, o->regions));
  }


  // check that our reinit function works correctly
#ifndef NDEBUG
  for (wall_index_t i = 0; i < p.get_wall_count(); i++) {
    Wall& w = p.get_wall(i);
    for (edge_index_t k = 0; k < EDGES_IN_TRIANGLE; k++) {
      Edge& e = w.edges[k];
      e.debug_check_values_are_uptodate(p);
    }
  }
#endif

  // --- back to object ---

  CHECK_PROPERTY(o->n_tiles == 0);
  CHECK_PROPERTY(o->n_occupied_tiles == 0);
  CHECK_PROPERTY(o->n_occupied_tiles == 0);
  CHECK_PROPERTY(t_matrix_to_mat4x4(o->t_matrix) == mat4x4(1) && "only identity matrix for now");
  // root->is_closed - not checked
  CHECK_PROPERTY(o->periodic_x == false);
  CHECK_PROPERTY(o->periodic_y == false);
  CHECK_PROPERTY(o->periodic_z == false);

  return true;
}


// cannot fail
void MCell3WorldConverter::create_diffusion_events() {
  assert(world->get_all_species().get_count() != 0 && "There must be at least 1 species");

  set<float_t> time_steps_set;
  for (auto &species : world->get_all_species().get_species_vector() ) {
    time_steps_set.insert(species.time_step);
  }

  for (float_t time_step : time_steps_set) {
    DiffuseReactEvent* event = new DiffuseReactEvent(world, time_step);
    event->event_time = TIME_SIMULATION_START;
    world->scheduler.schedule_event(event);
  }
}

static bool is_species_superclass(volume* s, species* spec) {
  return spec == s->all_mols || spec == s->all_volume_mols || spec == s->all_surface_mols;
}

bool MCell3WorldConverter::convert_species(volume* s) {

  // TODO_CONVERSION: many items are not checked
  for (int i = 0; i < s->n_species; i++) {
    species* spec = s->species_list[i];

    Species new_species;
    new_species.name = get_sym_name(spec->sym);
    new_species.D = spec->D;
    new_species.space_step = spec->space_step;
    new_species.time_step = spec->time_step;

    CHECK_PROPERTY((spec->flags & CANT_INITIATE) == 0); // RxnClass::compute_pb_factor assumes this

    // remove some flags for check that are known to work in all cases
    uint flags_check = spec->flags & ~REGION_PRESENT;

    if (!(
        is_species_superclass(s, spec)
        || flags_check == 0

        || flags_check == SPECIES_CPLX_MOL_FLAG_SURF
        || flags_check == SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE

        || flags_check == (SPECIES_FLAG_COUNT_ENCLOSED | COUNT_CONTENTS)

        || flags_check == SPECIES_FLAG_CAN_VOLVOL
        || flags_check == SPECIES_FLAG_CAN_VOLWALL
        || flags_check == SPECIES_FLAG_CAN_VOLSURF
        || flags_check == (SPECIES_FLAG_CAN_VOLVOL | SPECIES_FLAG_CAN_VOLSURF)
        || flags_check == (SPECIES_FLAG_CAN_VOLVOL | SPECIES_FLAG_CAN_VOLWALL)
        || flags_check == (SPECIES_FLAG_CAN_VOLWALL | SPECIES_FLAG_CAN_VOLSURF)

        || flags_check == (SPECIES_CPLX_MOL_FLAG_SURF | SPECIES_FLAG_CAN_SURFSURF)
        || flags_check == (SPECIES_CPLX_MOL_FLAG_SURF | SPECIES_FLAG_CAN_REGION_BORDER)

        || flags_check == (SPECIES_FLAG_CAN_VOLVOL | SPECIES_FLAG_CAN_VOLWALL | COUNT_CONTENTS | COUNT_ENCLOSED)
        || flags_check == (SPECIES_CPLX_MOL_FLAG_SURF | COUNT_CONTENTS)
        || flags_check == (SPECIES_CPLX_MOL_FLAG_SURF | SPECIES_FLAG_CAN_VOLVOL | COUNT_CONTENTS)
        || flags_check == (SPECIES_CPLX_MOL_FLAG_SURF | SPECIES_FLAG_CAN_SURFSURF | COUNT_CONTENTS)
        || flags_check == (SPECIES_CPLX_MOL_FLAG_SURF | SPECIES_FLAG_CAN_SURFSURF | SPECIES_FLAG_CAN_VOLWALL | COUNT_CONTENTS)
        || flags_check == (SPECIES_FLAG_CAN_VOLWALL | SPECIES_FLAG_COUNT_ENCLOSED | COUNT_CONTENTS)
        || flags_check == (SPECIES_CPLX_MOL_FLAG_SURF | SPECIES_FLAG_CAN_SURFSURF | CAN_SURFWALL | SPECIES_FLAG_CAN_REGION_BORDER)

      )) {
      mcell_log("Unsupported species flag for species %s: %s\n", new_species.name.c_str(), get_species_flags_string(spec->flags).c_str());
      CHECK_PROPERTY(false && "Flags listed in the message above are not supported yet");
    }

    // simply copy all the flags from mcell3
    new_species.set_flags(spec->flags);

    CHECK_PROPERTY(spec->n_deceased == 0);
    CHECK_PROPERTY(spec->cum_lifetime_seconds == 0);

    if ((spec->flags & IS_SURFACE) != 0) {
      if (spec->refl_mols != nullptr) {
        // just one type for now
        CHECK_PROPERTY(spec->refl_mols == nullptr || spec->refl_mols->next == nullptr);

        // reflective surface, seems that this information is transformed into reactions, so we do no need to store anything else
        CHECK_PROPERTY(spec->refl_mols->orient == ORIENTATION_NONE);
      }
    }
    else {
      CHECK_PROPERTY(spec->refl_mols == nullptr);
    }

    // CHECK_PROPERTY(spec->transp_mols == nullptr);
    // CHECK_PROPERTY(spec->absorb_mols == nullptr);
    CHECK_PROPERTY(spec->clamp_conc_mols == nullptr);

    // we must add a complex instance as the single molecule type in the new species
    // define a molecule type with no components
    MolType mol_type;
    mol_type.name = new_species.name; // name of the mol type is the same as for our species
    mol_type_id_t mol_type_id = world->bng_engine.get_data().find_or_add_molecule_type(mol_type);

    MolInstance mol_inst;
    mol_inst.mol_type_id = mol_type_id;
    if ((spec->flags & ON_GRID) == 0 && (spec->flags & IS_SURFACE) == 0) {  //FIXME: ALL_SURF_MOLS have wrong flag
      mol_inst.set_is_vol();
    }
    else if ((spec->flags & ON_GRID) != 0) {
      mol_inst.set_is_surf();
    }
    else if ((spec->flags & IS_SURFACE) != 0) {
      mol_inst.set_is_reactive_surface();
    }
    else {
      assert(false);
    }

    new_species.mol_instances.push_back(mol_inst);

    // and finally let's add our new species
    new_species.finalize();
    species_id_t new_species_id = world->get_all_species().find_or_add(new_species);

    // set all species 'superclasses' ids
    // these special species might be used in wall - surf|vol reactions
    if (spec == s->all_mols) {
      CHECK_PROPERTY(new_species.name == ALL_MOLECULES);
      world->get_all_species().set_all_molecules_species_id(new_species_id);
    }
    else if (spec == s->all_volume_mols) {
      CHECK_PROPERTY(new_species.name == ALL_VOLUME_MOLECULES);
      world->get_all_species().set_all_volume_molecules_species_id(new_species_id);
    }
    else if (spec == s->all_surface_mols) {
      CHECK_PROPERTY(new_species.name == ALL_SURFACE_MOLECULES);
      world->get_all_species().set_all_surface_molecules_species_id(new_species_id);
    }

    // map for other conversion steps
    mcell3_species_id_map[spec->species_id] = new_species_id;
  }

  return true;
}


// variable_rates_per_pathway is already resized to the number of pathways
static bool get_variable_rates(const t_func* prob_t, small_vector<small_vector<BNG::RxnRateInfo>>& variable_rates_per_pathway) {
  for (const t_func* curr = prob_t; curr != nullptr; curr = curr->next) {
    assert(curr->path >= 0 && curr->path < (int)variable_rates_per_pathway.size());

    BNG::RxnRateInfo info;
    info.time = curr->time;
    info.rate_constant = curr->value_from_file;

    assert(curr->path < (int)variable_rates_per_pathway.size());
    variable_rates_per_pathway[curr->path].push_back(info);
  }

  // make sure that they are sorted by time
  for (small_vector<BNG::RxnRateInfo>& rates: variable_rates_per_pathway) {
    sort(variable_rates_per_pathway.begin(), variable_rates_per_pathway.end());
  }

  return true;
}


bool MCell3WorldConverter::convert_single_reaction(const rxn *mcell3_rx) {
  // rx->next - handled in convert_reactions
  // rx->sym->name - ignored, name obtained from pathway

  //?? u_int n_reactants - obtained from pathways, might check it

  CHECK_PROPERTY(
      mcell3_rx->n_pathways >= 1
      || mcell3_rx->n_pathways == RX_REFLEC // reflections for surf mols
      || mcell3_rx->n_pathways == RX_TRANSP
  ); // limited for now

  assert(mcell3_rx->cum_probs != nullptr);

  // BNGTODO: create RXN classes!
  // rxn_class.max_fixed_p = mcell3_rx->max_fixed_p;
  // rxn_class.min_noreaction_p = mcell3_rx->min_noreaction_p;

  // ?? double pb_factor; /* Conversion factor from rxn rate to rxn probability (used for cooperativity) */

  // int *product_idx_aux - ignored, post-processing information
  // u_int *product_idx - ignored, post-processing information
  // struct species **players - ignored, might check it but will contain the same info as pathways

  // NFSIM struct species ***nfsim_players; /* a matrix of the nfsim elements associated with each path */
  // NFSIM short *geometries;         /* Geometries of reactants/products */
  // NFSIM short **nfsim_geometries;   /* geometries of the nfsim geometries associated with each path */

  CHECK_PROPERTY(mcell3_rx->n_occurred == 0);
  CHECK_PROPERTY(mcell3_rx->n_skipped == 0);

  small_vector<small_vector<BNG::RxnRateInfo>> variable_rates_per_pathway;
  if (mcell3_rx->n_pathways > 0) {
    variable_rates_per_pathway.resize(mcell3_rx->n_pathways);

    get_variable_rates(mcell3_rx->prob_t, variable_rates_per_pathway);
  }
  else {
    CHECK_PROPERTY(mcell3_rx->prob_t == nullptr);
  }

  // TODO_CONVERSION: pathway_info *info - magic_list, also some checks might be useful

  int pathway_index = 0;
  for (
      pathway* current_pathway = mcell3_rx->pathway_head;
      current_pathway != nullptr;
      current_pathway = current_pathway->next) {

    // --- pathway ---

    // there can be a single reaction for absorb, reflect and transparent reactions
    assert((pathway_index < mcell3_rx->n_pathways) || (pathway_index == 0 && mcell3_rx->n_pathways < 0));

    // -> pathway is renamed in MCell3 to reaction because pathway has a different meaning
    //    MCell3 reaction is reaction class
    RxnRule rxn;

    if (current_pathway->pathname != nullptr) {
      assert(current_pathway->pathname->sym != nullptr);
      rxn.name = get_sym_name(current_pathway->pathname->sym);
    }
    else if (mcell3_rx->sym != nullptr) {
      rxn.name = get_sym_name(mcell3_rx->sym);
    }
    else {
      rxn.name = NAME_NOT_SET;
    }

    if (mcell3_rx->prob_t == nullptr) {
      rxn.rate_constant = current_pathway->km;
    }
    else {
      CHECK_PROPERTY(mcell3_rx->n_pathways > 0);
      CHECK_PROPERTY(mcell3_rx->pb_factor != 0);

      // with variable rates, we may need to recompute the initial value for iteration 0
      float_t prob = (pathway_index == 0) ?
          mcell3_rx->cum_probs[0] : (mcell3_rx->cum_probs[pathway_index] - mcell3_rx->cum_probs[pathway_index-1]);

      rxn.rate_constant = prob / mcell3_rx->pb_factor;
    }

    if (!variable_rates_per_pathway.empty()) {
      assert(pathway_index < (int)variable_rates_per_pathway.size());
      rxn.variable_rates = variable_rates_per_pathway[pathway_index];
    }

    CHECK_PROPERTY(current_pathway->reactant1 != nullptr);
    CHECK_PROPERTY(current_pathway->orientation1 == 0 || current_pathway->orientation1 == 1 || current_pathway->orientation1 == -1);
    CHECK_PROPERTY(
        current_pathway->reactant2 == nullptr ||
        (current_pathway->orientation2 == 0 || current_pathway->orientation2 == 1 || current_pathway->orientation2 == -1)
    );
    CHECK_PROPERTY(
        current_pathway->reactant2 == nullptr ||
        (current_pathway->orientation3 == 0 || current_pathway->orientation3 == 1 || current_pathway->orientation3 == -1)
    );

    // reactants
    if (current_pathway->reactant1 != nullptr) {
      species_id_t reactant1_id = get_mcell4_species_id(current_pathway->reactant1->species_id);
      rxn.append_reactant(world->bng_engine.create_species_based_cplx_instance(reactant1_id, current_pathway->orientation1));

      if (current_pathway->reactant2 != nullptr) {
        species_id_t reactant2_id = get_mcell4_species_id(current_pathway->reactant2->species_id);
        rxn.append_reactant(world->bng_engine.create_species_based_cplx_instance(reactant2_id, current_pathway->orientation2));

        if (current_pathway->reactant3 != nullptr) {
          mcell_error("TODO_CONVERSION: reactions with 3 reactants are not supported");
        }
      }
      else {
        // reactant3 must be null if reactant2 is null
        assert(current_pathway->reactant3 == nullptr);
      }
    }
    else {
      assert(false && "No reactants?");
    }

    CHECK_PROPERTY(current_pathway->km_filename == nullptr);


    if (mcell3_rx->n_pathways == RX_ABSORB_REGION_BORDER) {
      CHECK_PROPERTY(current_pathway->flags == PATHW_ABSORP);
      rxn.type = RxnType::AbsorbRegionBorder;
    }
    else if (mcell3_rx->n_pathways == RX_TRANSP) {
      CHECK_PROPERTY(current_pathway->flags == PATHW_TRANSP);
      rxn.type = RxnType::Transparent;
    }
    else if (mcell3_rx->n_pathways == RX_REFLEC) {
      CHECK_PROPERTY(current_pathway->flags == PATHW_REFLEC);
      rxn.type = RxnType::Reflect;
    }
    else {
      CHECK_PROPERTY(current_pathway->flags == 0 || current_pathway->flags == PATHW_ABSORP);
      rxn.type = RxnType::Standard;

      // products
      product *product_ptr = current_pathway->product_head;
      while (product_ptr != nullptr) {
        CHECK_PROPERTY(product_ptr->orientation == 0 || product_ptr->orientation == 1 || product_ptr->orientation == -1);
        species_id_t product_id = get_mcell4_species_id(product_ptr->prod->species_id);
        rxn.append_product( world->bng_engine.create_species_based_cplx_instance(product_id, product_ptr->orientation) );

        product_ptr = product_ptr->next;
      }
    }

    pathway_index++;

    // add our reaction, reaction classes are created on-the-fly
    world->get_all_rxns().add_finalized_no_update(rxn);
  }



  return true;
}


bool MCell3WorldConverter::convert_rxns(volume* s) {

  rxn** reaction_hash = s->reaction_hash;
  int count = s->rx_hashsize;

  for (int i = 0; i < count; ++i) {
    rxn *rx = reaction_hash[i];
    while (rx != nullptr) {
      CHECK(convert_single_reaction(rx));
      rx = rx->next;
    }
  }

  return true;
}

// returns false if conversion failed
static RegionExprNode* create_release_region_terms_recursively(release_evaluator* expr, ReleaseEvent& event_data) {
  assert(expr != nullptr);
  RegionExprOperator new_op = RegionExprOperator::Invalid;

  uint op_masked = expr->op & REXP_MASK;

  switch (op_masked) {
    case REXP_UNION:
      new_op = RegionExprOperator::Union;
      break;
    case REXP_INTERSECTION:
      new_op = RegionExprOperator::Intersection;
      break;
    case REXP_SUBTRACTION:
      new_op = RegionExprOperator::Subtraction;
      break;

    case REXP_NO_OP:
      // do nothing, the mcell3 expression tree has no leaves and regions are stored directly
      break;

    default:
      assert(false && "Invalid region release_evaluator operator");
  }

  RegionExprNode* new_left;
  RegionExprNode* new_right;

  if ((expr->op & REXP_LEFT_REGION) != 0) {
    // left node is a region
    new_left = event_data.create_new_region_expr_node_leaf(get_sym_name(((region*)expr->left)->sym));
  }
  else if (new_op != RegionExprOperator::Invalid){
    // there is an operator so we assume that the left node is a subexpr
    new_left = create_release_region_terms_recursively((release_evaluator*)expr->left, event_data);
  }
  else {
    new_left = nullptr;
  }

  if ((expr->op & REXP_RIGHT_REGION) != 0) {
    new_right = event_data.create_new_region_expr_node_leaf(get_sym_name(((region*)expr->right)->sym));
  }
  else if (new_op != RegionExprOperator::Invalid) {
    new_right = create_release_region_terms_recursively((release_evaluator*)expr->right, event_data);
  }
  else {
    new_right = nullptr;
  }

  if (new_left != nullptr && new_right != nullptr) {
    assert(new_op != RegionExprOperator::Invalid);
    return event_data.create_new_region_expr_node_op(new_op, new_left, new_right);
  }
  else if (new_left != nullptr) {
    return new_left;
  }
  else if (new_right != nullptr) {
    return new_right;
  }
  else {
    assert(false);
    return nullptr;
  }
}

// FIXME: thiscould use a cleanup and maybe redesign
bool MCell3WorldConverter::convert_release_events(volume* s) {


  // -- schedule_helper -- (as volume.releaser)
  for (schedule_helper* releaser = s->releaser; releaser != NULL; releaser = releaser->next_scale) {

    //CHECK_PROPERTY(releaser->dt == 1); // it seems that wecan safely ignore dt
    //CHECK_PROPERTY(releaser->dt_1 == 1);
    // CHECK_PROPERTY(releaser->now == 0); // and now as well?
    //ok now: CHECK_PROPERTY(releaser->count == 1);
    CHECK_PROPERTY(releaser->index == 0);

    for (int i = -1; i < releaser->buf_len; i++) {
      for (abstract_element *aep = (i < 0) ? releaser->current : releaser->circ_buf_head[i]; aep != NULL; aep = aep->next) {

        ReleaseEvent event_data(world); // used only locally to capture the information

        // -- release_event_queue --
        release_event_queue *req = (release_event_queue *)aep;

        // mcell3 has several independent scheduling queues, this means that a release scheduled for
        // middle of of a timestep is executed at the beginning of a timestep, time is the number of iterations
        // therefore to get the iteration time, we can simply floor the value
        event_data.event_time = floor_to_multiple(req->event_time, 1);
        event_data.actual_release_time = req->event_time;

        // -- release_site --
        release_site_obj* rel_site = req->release_site;

        CHECK_PROPERTY(rel_site->release_shape == SHAPE_SPHERICAL || rel_site->release_shape == SHAPE_REGION);
        switch (rel_site->release_shape) {
          case SHAPE_SPHERICAL:
            event_data.release_shape = ReleaseShape::SPHERICAL;
            break;
          case SHAPE_REGION:
            event_data.release_shape = ReleaseShape::REGION;
            break;
          default:
            assert(false);
        }

        if (rel_site->region_data == nullptr) {
          assert(rel_site->location != nullptr);
          event_data.location = Vec3(*rel_site->location);
        }
        else if (rel_site->location == nullptr) {
          assert(rel_site->region_data != nullptr);
          release_region_data* region_data = rel_site->region_data;
          event_data.location = Vec3(POS_INVALID);

          // CHECK_PROPERTY(region_data->in_release == nullptr); // not sure what this means yet


          if (rel_site->region_data->cum_area_list != nullptr) {
            // surface molecules release onto region

            release_region_data* region_data = rel_site->region_data;
            for (int wall_i = 0; wall_i < region_data->n_walls_included; wall_i++) {

              wall* w = region_data->owners[region_data->obj_index[wall_i]]->wall_p[region_data->wall_index[wall_i]];

              PartitionWallIndexPair wall_index = get_mcell4_wall_index(w);

              event_data.cum_area_and_pwall_index_pairs.push_back(
                  CummAreaPWallIndexPair(region_data->cum_area_list[wall_i], wall_index)
              );
            }

            // also remember the expression, although the cum_area_and_pwall_index_pairs is what is currently used
            CHECK_PROPERTY(region_data->expression != nullptr);
            event_data.region_expr_root = create_release_region_terms_recursively(region_data->expression, event_data);
          }
          else {
            // volume or surface molecules release into region

            // this is quite limited for now, a single region is allowed
            CHECK_PROPERTY(region_data->cum_area_list == nullptr);
            CHECK_PROPERTY(region_data->wall_index == nullptr);
            CHECK_PROPERTY(region_data->n_objects == -1);
            CHECK_PROPERTY(region_data->owners == 0);
            // CHECK_PROPERTY(region_data->walls_per_obj == 0); not sure what this does

            event_data.region_llf = region_data->llf;
            event_data.region_urb = region_data->urb;

            CHECK_PROPERTY(region_data->expression != nullptr);
            event_data.region_expr_root = create_release_region_terms_recursively(region_data->expression, event_data);
          }
        }
        else {
          CHECK_PROPERTY(
              false
              && "So far supporting location for volume molecules and region for surface molecules"
          );
        }

        event_data.species_id = get_mcell4_species_id(rel_site->mol_type->species_id);
        assert(world->get_all_species().is_valid_id(event_data.species_id));

        CHECK_PROPERTY(rel_site->release_number_method == CONSTNUM || rel_site->release_number_method == DENSITYNUM);
        switch(rel_site->release_number_method) {
          case CONSTNUM:
            event_data.release_number_method = ReleaseNumberMethod::ConstNum;
            CHECK_PROPERTY(rel_site->concentration == 0);
            break;
          case DENSITYNUM:
            event_data.release_number_method = ReleaseNumberMethod::DensityNum;
            break;
          default:
            assert(false);
        }


        CHECK_PROPERTY(event_data.release_shape == ReleaseShape::REGION || rel_site->orientation == 0);
        event_data.orientation = rel_site->orientation;

        event_data.release_number = rel_site->release_number;

        CHECK_PROPERTY(rel_site->mean_diameter == 0); // temporary

        event_data.concentration = rel_site->concentration;

        CHECK_PROPERTY(rel_site->standard_deviation == 0); // temporary

        CHECK_PROPERTY(event_data.release_shape == ReleaseShape::REGION || rel_site->diameter != nullptr);
        if (rel_site->diameter != nullptr) {
          event_data.diameter = *rel_site->diameter;
        }
        else {
          event_data.diameter = Vec3(LENGTH_INVALID);
        }


        CHECK_PROPERTY(rel_site->mol_list == nullptr);
        CHECK_PROPERTY(rel_site->release_prob == 1); // temporary
        // rel_site->periodic_box - ignoring?

        event_data.release_site_name = rel_site->name;
        // rel_site->graph_pattern - ignored

        // -- release_event_queue -- (again)
        CHECK_PROPERTY(t_matrix_to_mat4x4(req->t_matrix) == mat4x4(1) && "only identity matrix for now");
        CHECK_PROPERTY(req->train_counter == 0);

        //release_pattern
        assert(rel_site->pattern != nullptr);
        release_pattern* rp = rel_site->pattern;
        if (rp->sym != nullptr) {
          event_data.release_pattern_name = get_sym_name(rp->sym);
        }
        else {
          event_data.release_pattern_name = NAME_NOT_SET;
        }

        assert(rp->delay == req->event_time && "Release pattern must specify the same delay as for which the event is scheduled");
        assert(rp->delay == req->train_high_time && "Initial train_high_time must be the same as delay");

        // schedule all the needed release events based on release pattern
        // note: there might be many of them but for now, we assume that not so many
        // maybe we will need to change it in a way so that the event schedules itself, but this was
        // a simpler solution for now
        float_t next_time = rp->delay;
        for (int train = 0; train < rp->number_of_trains; train++) {

          float_t train_start = rp->delay + train * rp->train_interval;
          // -EPS is needed to deal with precision issues even when we 
          // are strictly (<) comparing current and end time
          float_t train_end = train_start + rp->train_duration - EPS; 
          float_t current_time = train_start;
          while (current_time < train_end) {
            ReleaseEvent* event_to_schedule = new ReleaseEvent(world);
            *event_to_schedule = event_data;

            // similar as above for releases without release patterns, we need to set time to the beginning of the iteration
            event_to_schedule->event_time = floor_to_multiple(current_time, 1);
            event_to_schedule->actual_release_time =
                event_to_schedule->event_time + (event_data.actual_release_time - event_data.event_time);
            world->scheduler.schedule_event(event_to_schedule); // we always need to schedule a new instance

            current_time += rp->release_interval;
          }
        }
      }
    }

    // -- schedule_helper -- (again)
    CHECK_PROPERTY(releaser->current_count ==  0);
    CHECK_PROPERTY(releaser->current == nullptr);
    CHECK_PROPERTY(releaser->current_tail == nullptr);
    CHECK_PROPERTY(releaser->defunct_count == 0);
    CHECK_PROPERTY(releaser->error == 0);
    //CHECK_PROPERTY(releaser->depth == 0); // ignore

  }

  return true;
}


bool MCell3WorldConverter::convert_viz_output_events(volume* s) {

  // -- viz_output_block --
  viz_output_block* viz_blocks = s->viz_blocks;
  if (viz_blocks == nullptr) {
    return true; // no visualization data
  }
  CHECK_PROPERTY(viz_blocks->next == nullptr);
  CHECK_PROPERTY(viz_blocks->viz_mode == NO_VIZ_MODE || viz_blocks->viz_mode  == ASCII_MODE || viz_blocks->viz_mode == CELLBLENDER_MODE); // just checking valid values
  viz_mode_t viz_mode = viz_blocks->viz_mode;
  // CHECK_PROPERTY(viz_blocks->viz_output_flag == VIZ_ALL_MOLECULES); // ignored (for now?)

  uint_set<species_id_t> species_ids_to_visualize;
  assert((int)mcell3_species_id_map.size() == s->n_species);
  for (int i = 0; i < s->n_species; i++) {
    if (viz_blocks->species_viz_states[i] == INCLUDE_OBJ) {
      species_ids_to_visualize.insert_unique( get_mcell4_species_id(i) );
    }
    else if (viz_blocks->species_viz_states[i] == EXCLUDE_OBJ) {
      continue; // this species is not counted
    }
    else {
      assert(false);
    }
  }

  // CHECK_PROPERTY(viz_blocks->default_mol_state == INCLUDE_OBJ); // seems to be used only internally

  // -- frame_data_head --
  frame_data_list* frame_data_head = viz_blocks->frame_data_head;
  if (frame_data_head == nullptr) {
    // no events to log
    return true;
  }

  CHECK_PROPERTY(frame_data_head->next == nullptr);
  CHECK_PROPERTY(frame_data_head->list_type == OUTPUT_BY_ITERATION_LIST); // limited for now
  CHECK_PROPERTY(frame_data_head->type == ALL_MOL_DATA); // limited for now
  CHECK_PROPERTY(frame_data_head->viz_iteration == 0); // must be zero at sim. start

  // determine periodicity
  CHECK_PROPERTY(frame_data_head->n_viz_iterations > 0);

  float_t periodicity;
  float_t start_time;

  if (frame_data_head->n_viz_iterations > 1) {
    num_expr_list* iteration_ptr = frame_data_head->iteration_list;

    assert(iteration_ptr->next != nullptr);
    periodicity = iteration_ptr->next->value - iteration_ptr->value;

    // check that the first different in time is true for everything else
    num_expr_list* curr_viz_iteration_ptr = frame_data_head->curr_viz_iteration;
    for (long long i = 0; i < frame_data_head->n_viz_iterations; i++) {
      assert(iteration_ptr != nullptr);
      assert(curr_viz_iteration_ptr != nullptr);
      assert(iteration_ptr->value == curr_viz_iteration_ptr->value);

      if (iteration_ptr->next != nullptr) {
        CHECK_PROPERTY(cmp_eq(periodicity, iteration_ptr->next->value - iteration_ptr->value));
      }
      iteration_ptr = iteration_ptr->next;
      curr_viz_iteration_ptr = curr_viz_iteration_ptr->next;
    }
  }
  else {
    periodicity = 0;
  }

  VizOutputEvent* event = new VizOutputEvent(world);
  event->event_time = frame_data_head->iteration_list->value;
  event->periodicity_interval = periodicity;
  event->viz_mode = viz_mode;
  event->file_prefix_name = viz_blocks->file_prefix_name;
  event->species_ids_to_visualize = species_ids_to_visualize; // copying the whole map

  world->scheduler.schedule_event(event);

  return true;
}


static output_request* find_output_request_by_requester(volume* s, output_expression* expr) {
  for (
      output_request* req = s->output_request_head;
      req != nullptr;
      req = req->next) {

    // pointer comparison
    if (req->requester == expr) {
      return req;
    }
  }

  return nullptr;
}

// returns false if conversion failed
static bool find_output_requests_terms_recursively(
    volume* s, output_expression* expr, int sign, std::vector< pair<output_request*, int>>& requests_with_sign) {
  assert(sign == -1 || sign == 1);

  // operator must me one of +, -, #
  // # - cast output_request to expr

  switch (expr->oper) {
    case '#': {
      output_request* req = find_output_request_by_requester(s, expr);
      CHECK_PROPERTY(req != nullptr);
      requests_with_sign.push_back(make_pair(req, sign));
      }
      break;


    case '+':
    case '-': {
        bool left_ok = find_output_requests_terms_recursively(s, (output_expression*)expr->left, sign, requests_with_sign);
        if (!left_ok) {
          return false;
        }

        if (expr->oper == '-') {
          sign = -sign;
        }
        bool right_ok = find_output_requests_terms_recursively(s, (output_expression*)expr->right, sign, requests_with_sign);
        if (!right_ok) {
          return false;
        }
      }
      break;
    default:
      CHECK_PROPERTY(false && "Invalid output request operator");
  }

  return true;
}

static bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) {
      return false;
    }
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

bool MCell3WorldConverter::convert_mol_count_and_rxn_count_events(volume* s) {
  output_block* output_blocks = s->output_block_head;

  if (output_blocks == nullptr) {
    return true;
  }

  CHECK_PROPERTY(output_blocks->next == nullptr); // for now just one block

  MolOrRxnCountEvent* event = new MolOrRxnCountEvent(world);

  CHECK_PROPERTY(output_blocks->timer_type == OUTPUT_BY_STEP);

  CHECK_PROPERTY(output_blocks->t == 0);
  event->event_time = 0;
  event->periodicity_interval = output_blocks->step_time / world->config.time_unit;
  size_t buffer_size = output_blocks->buffersize;

  // we can check that the time_array contains expected values

  for (
      output_set* data_set = output_blocks->data_set_head;
      data_set != nullptr;
      data_set = data_set->next) {

    CHECK_PROPERTY(data_set->outfile_name != nullptr);
    // appending "4" to distinguish from
    count_buffer_id_t buffer_id = world->create_count_buffer(string("4") + data_set->outfile_name, buffer_size);

    // NOTE: FILE_SUBSTITUTE is interpreted in the same way as FILE_OVERWRITE
    CHECK_PROPERTY(data_set->file_flags == FILE_OVERWRITE || data_set->file_flags == FILE_SUBSTITUTE);
    CHECK_PROPERTY(data_set->chunk_count == 0);
    CHECK_PROPERTY(data_set->header_comment == nullptr);
    CHECK_PROPERTY(data_set->exact_time_flag == 1);

    output_column* column_head = data_set->column_head;
    CHECK_PROPERTY(column_head->next == nullptr);
    CHECK_PROPERTY(column_head->initial_value == 0);

    // TODO: we should better check column_head->expr contents

    // this information is split in a weird way in MCell3,
    // output request contains more information on what should be counted,
    // there is a way how to get to the output_set from the expression
    // but not how to the output_request
    // we simplify it to a list of terms with their sign
    // first is the output request, second is sign (-1 or +1)
    std::vector< pair<output_request*, int>> requests_with_sign;
    CHECK(find_output_requests_terms_recursively(s, column_head->expr, +1, requests_with_sign));

    MolOrRxnCountInfo info(buffer_id);

    for (pair<output_request*, int>& req_w_sign: requests_with_sign) {

      MolOrRxnCountTerm term;

      term.sign_in_expression = req_w_sign.second;

      output_request* req = req_w_sign.first;
      CHECK_PROPERTY(req != 0);

      int o = req->count_orientation;
      if (o == ORIENT_NOT_SET) {
        term.orientation = ORIENTATION_NOT_SET;
      }
      else {
        CHECK_PROPERTY(o == ORIENTATION_UP || o == ORIENTATION_NONE || o == ORIENTATION_DOWN);
        term.orientation = o;
      }

      // report type
      CHECK_PROPERTY(
          req->report_type == REPORT_CONTENTS ||
          req->report_type == REPORT_RXNS ||
          req->report_type == (REPORT_CONTENTS | REPORT_ENCLOSED) ||
          req->report_type == (REPORT_CONTENTS | REPORT_WORLD) ||
          req->report_type == (REPORT_RXNS | REPORT_WORLD) ||
          req->report_type == (REPORT_RXNS | REPORT_ENCLOSED)
      );

      bool count_mols_not_rxns;
      if ((req->report_type & REPORT_CONTENTS) != 0) {
        count_mols_not_rxns = true;
      }
      else if ((req->report_type & REPORT_RXNS) != 0) {
        count_mols_not_rxns = false;
      }
      else {
        assert(false);
        count_mols_not_rxns = true; // silence compilation warning
      }

      // count location
      // only whole geom object for now
      if ((req->report_type & REPORT_ENCLOSED) != 0) {

        term.type = count_mols_not_rxns ? CountType::EnclosedInObject : CountType::RxnCountInObject;
        
        string reg_name = get_sym_name(req->count_location);
        CHECK_PROPERTY(ends_with(reg_name, ",ALL")); // enclused in object must be whole objects

        const Region* reg = world->get_partition(0).find_region_by_name(reg_name);
        CHECK_PROPERTY(reg != nullptr);
        assert(reg->geometry_object_id != GEOMETRY_OBJECT_ID_INVALID);
        term.geometry_object_id = reg->geometry_object_id;

        GeometryObject& obj = world->get_partition(0).get_geometry_object(term.geometry_object_id);

        // set flag that we should include this object in counted volumes
        obj.is_counted_volume = true;
      }
      else {
        if (req->count_location == nullptr) {
          term.type = count_mols_not_rxns ? CountType::EnclosedInWorld : CountType::RxnCountInWorld;
        }
        else {
          CHECK_PROPERTY(req->report_type == REPORT_CONTENTS || req->report_type == REPORT_RXNS);

          string reg_name = get_sym_name(req->count_location);
          const Region* reg = world->get_partition(0).find_region_by_name(reg_name);
          term.region_id = reg->id;

          term.type = count_mols_not_rxns ? CountType::PresentOnSurfaceRegion : CountType::RxnCountOnSurfaceRegion;
        }
      }

      if (count_mols_not_rxns) {
        // count target (species)
        CHECK_PROPERTY(req->count_target != 0);
        string species_name = get_sym_name(req->count_target);

        term.species_id = world->get_all_species().find_by_name(species_name);
        CHECK_PROPERTY(term.species_id != SPECIES_ID_INVALID);

        // set a flag that these species are to be counted
        BNG::Species& species = world->get_all_species().get(term.species_id);
        species.set_flag(BNG::SPECIES_FLAG_COUNT_ENCLOSED);
      }
      else {
        // set that the reaction must be counted
        CHECK_PROPERTY(req->count_target != nullptr);
        string rxn_name = get_sym_name(req->count_target);

        BNG::RxnRule* rxn_rule = world->get_all_rxns().find_rxn_rule_by_name(rxn_name);
        CHECK_PROPERTY(rxn_rule != nullptr && "Rxn rule with name was not found");
        rxn_rule->set_is_counted();

        term.rxn_rule_id = rxn_rule->id;
      }

      info.terms.push_back(term);
    }
    event->add_mol_count_info(info);
  }

  world->scheduler.schedule_event(event);

  return true;
}

} // namespace mcell
