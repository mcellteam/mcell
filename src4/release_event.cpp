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

#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs.h"
#include "logging.h"

#include <iostream>

#include "defines.h"

#include "release_event.h"
#include "world.h"
#include "partition.h"
#include "datamodel_defines.h"

#include "collision_utils.inc"
#include "geometry_utils.inc"
#include "grid_utils.inc"

using namespace std;

namespace MCell {

string RegionExprNode::to_string(const bool for_datamodel) const {
  stringstream out;
  assert(op != RegionExprOperator::Invalid);

  if (op == RegionExprOperator::Leaf) {
    if (for_datamodel) {
      return DMUtil::get_object_w_region_name(region_name);
    }
    else {
      return region_name;
    }
  }

  assert(left != nullptr);
  out << "(";
  out << left->to_string(for_datamodel);

  switch(op) {
    case RegionExprOperator::Union:
      out << " + ";
      break;
    case RegionExprOperator::Intersect:
      out << " * ";
      break;
    case RegionExprOperator::Difference:
      out << " - ";
      break;
    default:
      assert(false);
  }
  out << right->to_string(for_datamodel);
  out << ")";
  return out.str();
}


void RegionExprNode::dump() const {
  cout << to_string();

}


ReleaseEvent::~ReleaseEvent() {
  for (RegionExprNode* expr_node: all_region_expr_nodes) {
    delete expr_node;
  }
}


RegionExprNode* ReleaseEvent::create_new_region_expr_node_leaf(const std::string region_name) {
  RegionExprNode* res = new RegionExprNode;
  res->op = RegionExprOperator::Leaf;
  res->region_name = region_name;
  return res;
}


RegionExprNode* ReleaseEvent::create_new_region_expr_node_op(
    const RegionExprOperator op, RegionExprNode* left, RegionExprNode* right) {

  RegionExprNode* res = new RegionExprNode;
  res->op = op;
  res->left = left;
  res->right = right;
  return res;
}


void ReleaseEvent::dump(const string ind) const {
  cout << "Release event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "name: \t\t" << release_site_name << " [string] \t\t\n";
  cout << ind2 << "species_id: \t\t" << species_id << " [species_id_t] \t\t\n";
  cout << ind2 << "actual_release_time: \t\t" << location << " [float_t] \t\t\n";
  cout << ind2 << "release_number_method: \t\t" << (int)release_number_method << " [ReleaseNumberMethod] \t\t\n";
  cout << ind2 << "release_number: \t\t" << release_number << " [uint] \t\t\n";
  cout << ind2 << "concentration: \t\t" << concentration << " [float_t] \t\t\n";
  cout << ind2 << "orientation: \t\t" << orientation << " [float_t] \t\t\n";
  cout << ind2 << "release_shape: \t\t" << (int)release_shape << " [ReleaseShape] \t\t\n";
  cout << ind2 << "location: \t\t" << location << " [Vec3] \t\t\n";
  cout << ind2 << "diameter: \t\t" << diameter << " [Vec3] \t\t\n";

  // TODO:
  //cum_area_and_pwall_index_pairs
  //all_region_expr_nodes

  cout << ind2 << "region_llf: \t\t" << region_llf << " [Vec3] \t\t\n";
  cout << ind2 << "region_urb: \t\t" << region_urb << " [Vec3] \t\t\n";

  if (region_expr_root != nullptr) {
    region_expr_root->dump();
  }
}


void ReleaseEvent::to_data_model(Json::Value& mcell_node) const {

  if (event_time != 0 || actual_release_time != 0) {
    CONVERSION_UNSUPPORTED("Release event " + release_site_name + " starts at time different from 0, this is not supported yet.");
  }

  // these items might already exist
  Json::Value& release_sites = mcell_node[KEY_RELEASE_SITES];
  DMUtil::json_add_version(release_sites, JSON_DM_VERSION_1638);
  Json::Value& release_site_list = release_sites[KEY_RELEASE_SITE_LIST];

  Json::Value release_site;
  DMUtil::json_add_version(release_site, JSON_DM_VERSION_1330);
  release_site[KEY_DESCRIPTION] = "";
  release_site[KEY_POINTS_LIST] = Json::Value(Json::arrayValue); // not sure, empty array
  release_site[KEY_NAME] = DMUtil::remove_obj_name_prefix(release_site_name);
  release_site[KEY_MOLECULE] = world->get_all_species().get(species_id).name;
  release_site[KEY_ORIENT] = DMUtil::orientation_to_str(orientation);

  // how many to release
  switch (release_number_method) {
    case ReleaseNumberMethod::ConstNum:
      release_site[KEY_QUANTITY_TYPE] = VALUE_NUMBER_TO_RELEASE;
      release_site[KEY_QUANTITY] = to_string(release_number);
      break;
    case ReleaseNumberMethod::GaussNum:
      release_site[KEY_QUANTITY_TYPE] = VALUE_GAUSSIAN_RELEASE_NUMBER;
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has unsupported release_number_method GaussNum.");
      break;
    case ReleaseNumberMethod::VolNum:
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has unsupported release_number_method VolNum.");
      break;
    case ReleaseNumberMethod::ConcNum:
    case ReleaseNumberMethod::DensityNum:
      release_site[KEY_QUANTITY_TYPE] = VALUE_DENSITY;
      release_site[KEY_QUANTITY] = to_string(concentration);
      break;
    default:
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has invalid release_number_method.");
      break;
  }

  release_site[KEY_PATTERN] = "";
  release_site[KEY_STDDEV] = "0"; // TODO
  release_site[KEY_RELEASE_PROBABILITY] = DMUtil::f_to_string(1.0);  // only 1 for now

  // where to release
  switch (release_shape) {
    case ReleaseShape::SPHERICAL:
      release_site[KEY_SHAPE] = VALUE_SPHERICAL;
      break;
    case ReleaseShape::SPHERICAL_SHELL:
      release_site[KEY_SHAPE] = VALUE_SPHERICAL_SHELL;
      break;
    case ReleaseShape::REGION:
      release_site[KEY_SHAPE] = "OBJECT";
      release_site[KEY_OBJECT_EXPR] = region_expr_root->to_string(true);
      break;
    default:
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has shape different from SHPERE and OBJECT.");
      break;
  }

  if (release_shape != ReleaseShape::REGION) {
    release_site[KEY_LOCATION_X] = DMUtil::f_to_string(location.x * world->config.length_unit);
    release_site[KEY_LOCATION_Y] = DMUtil::f_to_string(location.y * world->config.length_unit);
    release_site[KEY_LOCATION_Z] = DMUtil::f_to_string(location.z * world->config.length_unit);

    CONVERSION_CHECK(diameter.x == diameter.y && diameter.y == diameter.z, "Not sure if datamodel supports different diameters.");
    release_site[KEY_SITE_DIAMETER] = diameter.x * world->config.length_unit;
  }

  release_site_list.append(release_site);
}


static void check_max_release_count(double num_to_release, const std::string& name) {
  int num = (int)num_to_release;
  if (num < 0 || num > INT_MAX) {
    mcell_error(
        "Release site '%s' tries to release more than INT_MAX (2147483647) molecules.",
        name.c_str());
  }
}


uint ReleaseEvent::calculate_number_to_release() {

  switch (release_number_method) {
    case ReleaseNumberMethod::ConstNum:
      return release_number;

    case ReleaseNumberMethod::ConcNum:
      if (diameter == Vec3(LENGTH_INVALID)) {
        return 0;
      }
      else {
        float_t vol;
        switch (release_shape) {
        case ReleaseShape::SPHERICAL:
        //case ReleaseShape::ELLIPTIC:
          vol = (1.0 / 6.0) * MY_PI * diameter.x * diameter.y * diameter.z;
          break;
        /*case SHAPE_RECTANGULAR:
        case SHAPE_CUBIC:
          vol = rso->diameter->x * rso->diameter->y * rso->diameter->z;
          break;*/

        case ReleaseShape::SPHERICAL_SHELL:
          mcell_error("Release site \"%s\" tries to release a concentration on a "
                      "spherical shell.", release_site_name.c_str());
          break;

        default:
          mcell_internal_error("Release by concentration on invalid release site "
                               "shape (%d) for release site \"%s\".",
                               (int)release_shape, release_site_name.c_str());
          break;
        }
        assert(concentration != FLT_INVALID);
        float_t num_to_release =
            N_AV * 1e-15 * concentration * vol * pow_f(world->config.length_unit, 3) + 0.5;
        check_max_release_count(num_to_release, release_site_name);
        return (uint)num_to_release;
      }
      break;

    case ReleaseNumberMethod::DensityNum: {
        // computed in release_onto_regions in MCell3
        assert(!cum_area_and_pwall_index_pairs.empty());
        float_t max_A = cum_area_and_pwall_index_pairs.back().first;
        float_t est_sites_avail = (int)max_A;

        assert(concentration != FLT_INVALID);
        float_t num_to_release = concentration * est_sites_avail / world->config.grid_density;
        check_max_release_count(num_to_release, release_site_name);
        return (uint)num_to_release;
      }
      break;

    default:
      assert(false);
      return 0;
  }
}


// NOTE: maybe a template will be needed for this function, used a lot in mcell3
static size_t cum_area_bisect_high(const vector<CummAreaPWallIndexPair>& array, float_t val) {
  size_t low = 0;
  size_t hi = array.size() - 1;
  size_t mid = 0;

  while (hi - low > 1) {
    mid = (hi + low) / 2;
    if (array[mid].first > val) {
      hi = mid;
    } else {
      low = mid;
    }
  }

  if (array[low].first > val)
  {
    return low;
  }
  else {
    return hi;
  }
}


void ReleaseEvent::place_single_molecule_onto_grid(Partition& p, Wall& wall, tile_index_t tile_index) {

  Vec2 pos_on_wall;
  if (p.config.randomize_smol_pos) {
    pos_on_wall = GridUtil::grid2uv_random(wall, tile_index, world->rng);
  }
  else {
    pos_on_wall = GridUtil::grid2uv(wall, tile_index);
  }

  Molecule& new_sm = p.add_surface_molecule(
      Molecule(MOLECULE_ID_INVALID, species_id, pos_on_wall, get_release_delay_time())
  );

  new_sm.s.wall_index = wall.index;
  new_sm.s.orientation = orientation;

  new_sm.s.grid_tile_index = tile_index;
  wall.grid.set_molecule_tile(tile_index, new_sm.id);

  new_sm.flags = ACT_DIFFUSE | IN_SURFACE;
  new_sm.set_flag(MOLECULE_FLAG_SURF);
  new_sm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

#ifdef DEBUG_RELEASES
  new_sm.dump(p, "Released vm:", "", p.stats.get_current_iteration(), actual_release_time, true);
#endif
}


void ReleaseEvent::release_onto_regions(uint computed_release_number) {
  int success = 0, failure = 0;
  float_t seek_cost = 0;

  // TODO_LATER: for now we are assuming that we have just a single partition
  // and releases do not cross partition boundary

  assert(!cum_area_and_pwall_index_pairs.empty());
  float_t total_area = cum_area_and_pwall_index_pairs.back().first;
  float_t est_sites_avail = (int)total_area;
  const float_t rel_list_gen_cost = 10.0; /* Just a guess */
  float_t pick_cost = rel_list_gen_cost * est_sites_avail;

  uint n = computed_release_number;

  const int too_many_failures = 10; /* Just a guess */
  while (n > 0) {
    if (failure >= success + too_many_failures) {
      seek_cost =
          n * (((double)(success + failure + 2)) / ((double)(success + 1)));
    }

    if (seek_cost < pick_cost) {
      float_t A = rng_dbl(&world->rng) * total_area;
      size_t cum_area_index = cum_area_bisect_high(cum_area_and_pwall_index_pairs, A);
      PartitionWallIndexPair pw = cum_area_and_pwall_index_pairs[cum_area_index].second;
      Partition& p = world->get_partition(pw.first);
      Wall& wall = p.get_wall(pw.second);

      if (!wall.has_initialized_grid()) {
        wall.initialize_grid(p); // sets wall's grid_index
      }

      Grid& grid = wall.grid;

      // get the random number for the current wall
      if (cum_area_index != 0) {
        A -= cum_area_and_pwall_index_pairs[cum_area_index - 1].first;
      }

      tile_index_t tile_index = (grid.num_tiles_along_axis * grid.num_tiles_along_axis) * (A / wall.area);
      if (tile_index >= grid.num_tiles) {
        tile_index = grid.num_tiles - 1;
      }

      if (grid.get_molecule_on_tile(tile_index) != MOLECULE_ID_INVALID) {
        failure++;
        continue;
      }

      place_single_molecule_onto_grid(p, wall, tile_index);

      success++;
      n--;
    }
    else {
      assert(false && "Recovery from too many failures during surf mol release is not implemented yet");
    }
  }
}


static bool is_point_inside_region_expr_recursively(const Partition& p, const Vec3& pos, const RegionExprNode* region_expr_node) {
  assert(region_expr_node->op != RegionExprOperator::Invalid);

  if (region_expr_node->op == RegionExprOperator::Leaf) {
    const Region* reg = p.find_region_by_name(region_expr_node->region_name);
    assert(reg != nullptr && "Region for release must exist");
    return CollisionUtil::is_point_inside_region(p, pos, *reg);
  }

  bool satisfies_l = is_point_inside_region_expr_recursively(p, pos, region_expr_node->left);
  bool satisfies_r = is_point_inside_region_expr_recursively(p, pos, region_expr_node->right);

  switch (region_expr_node->op) {
    case RegionExprOperator::Union:
      return satisfies_l || satisfies_r;
    case RegionExprOperator::Intersect:
      return satisfies_l && satisfies_r;
    case RegionExprOperator::Difference:
      return satisfies_l && !satisfies_r;
    default:
      assert(false);
      return false;
  }
}


void ReleaseEvent::release_inside_regions(uint computed_release_number) {

  assert(region_expr_root != nullptr);

  Partition& p = world->get_partition(0);

  /*if (rso->release_number_method == CCNNUM) {
    n = num_vol_mols_from_conc(rso, state->length_unit, &exactNumber);
  }*/

  int n = computed_release_number;

  while (n > 0) {
    Vec3 pos;
    pos.x = region_llf.x + (region_urb.x - region_llf.x) * rng_dbl(&world->rng);
    pos.y = region_llf.y + (region_urb.y - region_llf.y) * rng_dbl(&world->rng);
    pos.z = region_llf.z + (region_urb.z - region_llf.z) * rng_dbl(&world->rng);

    if (!is_point_inside_region_expr_recursively(p, pos, region_expr_root)) {
      /*if (rso->release_number_method == CCNNUM && !exactNumber)
        n--;*/
      continue;
    }

    // TODO_LATER: location can be close to a partition boundary, we might need to release to a different partition
    Molecule& new_vm = p.add_volume_molecule(
        Molecule(MOLECULE_ID_INVALID, species_id, pos, get_release_delay_time())
    );
    new_vm.flags = IN_VOLUME | ACT_DIFFUSE;
    new_vm.set_flag(MOLECULE_FLAG_VOL);
    new_vm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

    n--;

#ifdef DEBUG_RELEASES
    new_vm.dump(p, "Released vm:", "", p.stats.get_current_iteration(), actual_release_time, true);
#endif
  }
}


void ReleaseEvent::release_ellipsoid_or_rectcuboid(uint computed_release_number) {

  Partition& p = world->get_partition(world->get_or_add_partition_index(location));
  float_t time_step = world->get_all_species().get(species_id).time_step;

  const int is_spheroidal = (release_shape == ReleaseShape::SPHERICAL ||
                             /*release_shape == SHAPE_ELLIPTIC ||*/
                             release_shape == ReleaseShape::SPHERICAL_SHELL);

  for (uint i = 0; i < computed_release_number; i++) {
    Vec3 pos;
    do /* Pick values in unit square, toss if not in unit circle */
    {
      pos.x = (rng_dbl(&world->rng) - 0.5);
      pos.y = (rng_dbl(&world->rng) - 0.5);
      pos.z = (rng_dbl(&world->rng) - 0.5);
    } while (is_spheroidal && len3_squared(pos) >= 0.25);

    if (release_shape == ReleaseShape::SPHERICAL_SHELL) {
      float_t r = sqrt(len3_squared(pos)) * 2.0;
      if (r == 0.0) {
        pos = Vec3(0.0, 0.0, 0.5);
      } else {
        pos /= r;
      }
    }

    float_t base_location[1][4];
    base_location[0][0] = pos.x * diameter.x + location.x;
    base_location[0][1] = pos.y * diameter.y + location.y;
    base_location[0][2] = pos.z * diameter.z + location.z;
    base_location[0][3] = 1;

    // TODO_LATER: t_matrix can be only identity matrix for now, also use glm matrix mult.
    // mult_matrix(location, req->t_matrix, location, 1, 4, 4);

    Vec3 molecule_location;
    molecule_location.x = base_location[0][0];
    molecule_location.y = base_location[0][1];
    molecule_location.z = base_location[0][2];

    // TODO_LATER: location can be close to a partition boundary, we might need to release to a different partition
    Molecule& new_vm = p.add_volume_molecule(
        Molecule(MOLECULE_ID_INVALID, species_id, molecule_location, get_release_delay_time())
    );
    new_vm.flags = IN_VOLUME | ACT_DIFFUSE;
    new_vm.set_flag(MOLECULE_FLAG_VOL);
    new_vm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);
#ifdef DEBUG_RELEASES
    new_vm.dump(p, "Released vm:", "", p.stats.get_current_iteration(), actual_release_time, true);
#endif

  }
}


void ReleaseEvent::step() {
  // for now, let's simply release 'release_number' of molecules of 'species_id'
  // at 'location'

  uint number = calculate_number_to_release();

  const BNG::Species& species = world->get_all_species().get(species_id);

  if (release_shape == ReleaseShape::REGION) {
    if (species.is_surf()) {
      release_onto_regions(number);
    }
    else {
      release_inside_regions(number);
    }
  }
  else if (release_shape == ReleaseShape::SPHERICAL) {
    assert(diameter.is_valid());
    release_ellipsoid_or_rectcuboid(number);
  }
  else {
    assert(false);
  }

  cout
    << "Released " << number << " " << species.name << " from \"" << release_site_name << "\""
    << " at iteration " << world->get_current_iteration() << ".\n";
}


} // namespace mcell


