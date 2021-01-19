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
#include "mcell_structs_shared.h"
#include "logging.h"

#include <iostream>

#include "defines.h"

#include "release_event.h"
#include "world.h"
#include "partition.h"
#include "diffuse_react_event.h"
#include "datamodel_defines.h"

#include "geometry_utils.h"
#include "geometry_utils.inc"
#include "collision_utils.inc"
#include "grid_utils.inc"


using namespace std;

namespace MCell {

void dump_cumm_area_and_pwall_index_pairs(
    const std::vector<CummAreaPWallIndexPair>& cumm_area_and_pwall_index_pairs, const std::string ind) {
  cout << ind << "cumm_area_and_pwall_index_pairs:\n";
  size_t max = cumm_area_and_pwall_index_pairs.size();
#ifdef NDEBUG
  if (max > 4) {
    cout << ind << "Printing only the first 4 of " << max << "\n";
    max = 4;
  }
#endif

  for (size_t i = 0; i < max; i++) {
    cout << ind << i << ": ";
    const CummAreaPWallIndexPair& area_wall = cumm_area_and_pwall_index_pairs[i];
    cout <<
        "area: " << area_wall.first << ", partition_id: " << area_wall.second.first <<
        ", wall_index: " << area_wall.second.first << "\n";
  }
}


size_t cum_area_bisect_high(const std::vector<CummAreaPWallIndexPair>& array, float_t val) {
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


string RegionExprNode::to_string(const World* world, const bool for_datamodel) const {
  stringstream out;
  assert(op != RegionExprOperator::Invalid);

  if (op == RegionExprOperator::Leaf) {
    const string& region_name = world->get_region(region_id).name;
    if (for_datamodel) {
      return DMUtil::get_object_w_region_name(region_name);
    }
    else {
      return region_name;
    }
  }

  assert(left != nullptr);
  out << "(";
  out << left->to_string(world, for_datamodel);

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
  out << right->to_string(world, for_datamodel);
  out << ")";
  return out.str();
}


void RegionExprNode::dump(const World* world) const {
  cout << to_string(world);

}


ReleaseEvent::~ReleaseEvent() {
  for (RegionExprNode* expr_node: all_region_expr_nodes) {
    delete expr_node;
  }
}


bool ReleaseEvent::update_event_time_for_next_scheduled_time() {
  // see https://mcell.org/tutorials/_images/plot.png

  assert(number_of_trains == NUMBER_OF_TRAINS_UNLIMITED || current_train_from_0 <= number_of_trains);

  // did we process all trains?
  if (current_train_from_0 == number_of_trains) {
    return false;
  }

  // set the new times for schedule
  actual_release_time =
      delay +
      current_train_from_0 * train_interval +
      current_release_in_train_from_0 * release_interval;

  // floor_to_multiple adds EPS to the value before flooring
  event_time = floor_to_multiple(actual_release_time, 1);

  // increment release counter
  current_release_in_train_from_0++;

  // should we start a new train next time this method is called?
  int number_of_releases_per_train = get_num_releases_per_train();
  assert(number_of_releases_per_train >= 1);
  assert(current_release_in_train_from_0 <= number_of_releases_per_train);

  if (current_release_in_train_from_0 == number_of_releases_per_train) {
    current_train_from_0++;
    current_release_in_train_from_0 = 0;
  }

  return true;
}


RegionExprNode* ReleaseEvent::create_new_region_expr_node_leaf(const region_id_t region_id) {
  RegionExprNode* res = new RegionExprNode;
  res->op = RegionExprOperator::Leaf;
  res->region_id = region_id;
  all_region_expr_nodes.push_back(res);
  return res;
}


RegionExprNode* ReleaseEvent::create_new_region_expr_node_op(
    const RegionExprOperator op, RegionExprNode* left, RegionExprNode* right) {

  RegionExprNode* res = new RegionExprNode;
  res->op = op;
  res->left = left;
  res->right = right;
  all_region_expr_nodes.push_back(res);
  return res;
}

#define CASE_TO_STR(name) case name: return #name

static const char* release_shape_to_str(const ReleaseShape s) {
  switch (s) {
    CASE_TO_STR(ReleaseShape::UNDEFINED);
    CASE_TO_STR(ReleaseShape::SPHERICAL);
    CASE_TO_STR(ReleaseShape::SPHERICAL_SHELL);
    CASE_TO_STR(ReleaseShape::REGION);
    CASE_TO_STR(ReleaseShape::LIST);
    CASE_TO_STR(ReleaseShape::INITIAL_SURF_REGION);
    default: return "invalid_release_shape";
  }
}


static const char* release_number_method_to_str(const ReleaseNumberMethod m) {
  switch (m) {
    CASE_TO_STR(ReleaseNumberMethod::Invalid);
    CASE_TO_STR(ReleaseNumberMethod::ConstNum);
    CASE_TO_STR(ReleaseNumberMethod::GaussNum);
    CASE_TO_STR(ReleaseNumberMethod::VolNum);
    CASE_TO_STR(ReleaseNumberMethod::ConcentrationNum);
    CASE_TO_STR(ReleaseNumberMethod::DensityNum);
    default: return "invalid_release_number_method";
  }
}


void ReleaseEvent::dump(const string ind) const {
  cout << "Release event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "name: \t\t" << release_site_name << " [string]\n";
  cout << ind2 << "species_id: \t\t" << species_id << " [species_id_t]\n";
  cout << ind2 << "actual_release_time: \t\t" << actual_release_time << " [float_t]\n";
  cout << ind2 << "release_number_method: \t\t" << release_number_method_to_str(release_number_method) << " [ReleaseNumberMethod]\n";
  cout << ind2 << "release_number: \t\t" << release_number << " [uint]\n";
  cout << ind2 << "concentration: \t\t" << concentration << " [float_t]\n";
  cout << ind2 << "orientation: \t\t" << orientation << " [float_t]\n";
  cout << ind2 << "release_shape: \t\t" << release_shape_to_str(release_shape) << " [ReleaseShape]\n";
  cout << ind2 << "location: \t\t" << location << " [Vec3]\n";
  cout << ind2 << "diameter: \t\t" << diameter << " [Vec3]\n";

  dump_cumm_area_and_pwall_index_pairs(cumm_area_and_pwall_index_pairs, ind2);

  cout << ind2 << "region_llf: \t\t" << region_llf << " [Vec3]\n";
  cout << ind2 << "region_urb: \t\t" << region_urb << " [Vec3]\n";

  if (region_expr_root != nullptr) {
    region_expr_root->dump(world);
  }

  cout << ind2 << "delay: \t\t" << delay << " [float_t]\n";
  cout << ind2 << "number_of_trains: \t\t" << number_of_trains << " [uint]\n";
  cout << ind2 << "train_interval: \t\t" << train_interval << " [float_t]\n";
  cout << ind2 << "train_duration: \t\t" << train_duration << " [float_t]\n";
  cout << ind2 << "release_interval: \t\t" << release_interval << " [float_t]\n";
}


// returns the name of the release pattern,
// empty string if release pattern is not needed
std::string ReleaseEvent::release_pattern_to_data_model(Json::Value& mcell_node) const {
  if (!needs_release_pattern()) {
    // no release pattern is needed
    return "";
  }

  Json::Value& define_release_patterns = mcell_node[KEY_DEFINE_RELEASE_PATTERNS];
  // version might be already there
  if (define_release_patterns.isMember(KEY_DATA_MODEL_VERSION)) {
    DMUtil::add_version(define_release_patterns, VER_DM_2014_10_24_1638);
  }
  Json::Value& release_pattern_list = define_release_patterns[KEY_RELEASE_PATTERN_LIST];

  string name;
  if (release_pattern_name != "") {
    // should be usually set when rel pat is needed
    name = release_pattern_name;
  }
  else {
    name = RELEASE_PATTERN_PREFIX + release_site_name;
  }

  // the current pattern may be already present, skip it in this case,
  // release patterns must have unique names, checked while MDL is parsed
  // TODO: is this true also for MCell4?
  for (Json::ArrayIndex i = 0; i < release_pattern_list.size(); i++) {
    if (release_pattern_list[i][KEY_NAME].asString() == release_pattern_name) {
      return name;
    }
  }

  Json::Value release_pattern_item;

  DMUtil::add_version(release_pattern_item, VER_DM_2018_01_11_1330);

  release_pattern_item[KEY_DELAY] = f_to_str(delay * world->config.time_unit);
  release_pattern_item[KEY_NUMBER_OF_TRAINS] = to_string(number_of_trains);
  release_pattern_item[KEY_TRAIN_INTERVAL] = f_to_str(train_interval * world->config.time_unit);
  release_pattern_item[KEY_TRAIN_DURATION] = f_to_str(train_duration * world->config.time_unit);
  release_pattern_item[KEY_RELEASE_INTERVAL] = f_to_str(release_interval * world->config.time_unit);


  release_pattern_item[KEY_NAME] = name;
  release_pattern_item[KEY_DESCRIPTION] = "";

  release_pattern_list.append(release_pattern_item);

  return name;
}


void ReleaseEvent::to_data_model_as_one_release_site(
    Json::Value& mcell_node,
    const species_id_t species_id_override,
    const orientation_t orientation_override,
    // points_list indices are set only when
    // release_shape == ReleaseShape::LIST
    const string& name_override,
    const uint points_list_begin_index,
    const uint points_list_end_index
) const {

  // these items might already exist
  Json::Value& release_sites = mcell_node[KEY_RELEASE_SITES];
  DMUtil::add_version(release_sites, VER_DM_2014_10_24_1638);
  Json::Value& release_site_list = release_sites[KEY_RELEASE_SITE_LIST];

  Json::Value release_site;
  DMUtil::add_version(release_site, VER_DM_2018_01_11_1330);
  release_site[KEY_DESCRIPTION] = "";
  release_site[KEY_NAME] = DMUtil::remove_obj_name_prefix(name_override);
  release_site[KEY_MOLECULE] = world->get_all_species().get(species_id_override).name;
  release_site[KEY_ORIENT] = DMUtil::orientation_to_str(orientation_override);

  // how many to release
  switch (release_number_method) {
    case ReleaseNumberMethod::ConstNum:
      release_site[KEY_QUANTITY_TYPE] = VALUE_NUMBER_TO_RELEASE;
      if (release_shape != ReleaseShape::LIST) {
        release_site[KEY_QUANTITY] = to_string(release_number);
      }
      else {
        release_site[KEY_QUANTITY] = ""; // number for release of LIST is given by the number of points
      }
      break;
    case ReleaseNumberMethod::GaussNum:
      release_site[KEY_QUANTITY_TYPE] = VALUE_GAUSSIAN_RELEASE_NUMBER;
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has unsupported release_number_method GaussNum.");
      break;
    case ReleaseNumberMethod::VolNum:
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has unsupported release_number_method VolNum.");
      break;
    case ReleaseNumberMethod::ConcentrationNum:
    case ReleaseNumberMethod::DensityNum:
      release_site[KEY_QUANTITY_TYPE] = VALUE_DENSITY;
      release_site[KEY_QUANTITY] = f_to_str(concentration);
      break;
    default:
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has invalid release_number_method.");
      break;
  }


  string data_model_release_pattern_name = release_pattern_to_data_model(mcell_node);
  release_site[KEY_PATTERN] = data_model_release_pattern_name;

  release_site[KEY_STDDEV] = "0"; // TODO
  release_site[KEY_RELEASE_PROBABILITY] = f_to_str(1.0);  // only 1 for now

  // where to release
  switch (release_shape) {
    case ReleaseShape::SPHERICAL:
      release_site[KEY_SHAPE] = VALUE_SPHERICAL;
      break;
    case ReleaseShape::SPHERICAL_SHELL:
      release_site[KEY_SHAPE] = VALUE_SPHERICAL_SHELL;
      break;
    case ReleaseShape::REGION:
      release_site[KEY_SHAPE] = VALUE_OBJECT;
      release_site[KEY_OBJECT_EXPR] = region_expr_root->to_string(world, true);
      break;
    case ReleaseShape::LIST:
      release_site[KEY_SHAPE] = VALUE_LIST;
      break;
    default:
      CONVERSION_UNSUPPORTED("Release event " + release_site_name + " has shape different from SHPERE and OBJECT.");
      break;
  }

  if (release_shape != ReleaseShape::REGION) {
    if (release_shape != ReleaseShape::LIST) {
      release_site[KEY_LOCATION_X] = f_to_str(location.x * world->config.length_unit);
      release_site[KEY_LOCATION_Y] = f_to_str(location.y * world->config.length_unit);
      release_site[KEY_LOCATION_Z] = f_to_str(location.z * world->config.length_unit);
    }
    else {
      assert(points_list_begin_index < points_list_end_index);
      Json::Value& points_list = release_site[KEY_POINTS_LIST];
      for (uint i = points_list_begin_index; i < points_list_end_index; i++) {
        assert(i < molecule_list.size());
        Vec3 pos_scaled = molecule_list[i].pos * Vec3(world->config.length_unit);
        DMUtil::append_triplet(points_list, pos_scaled.x, pos_scaled.y, pos_scaled.z);
      }
    }

    CONVERSION_CHECK(diameter.x == diameter.y && diameter.y == diameter.z, "Datamodel does not support different diameters.");
    release_site[KEY_SITE_DIAMETER] = f_to_str(diameter.x * world->config.length_unit);
  }

  release_site_list.append(release_site);
}


void ReleaseEvent::to_data_model(Json::Value& mcell_node) const {

  if (event_time != 0 && !needs_release_pattern()) {
    // the MCell4 API supports this, but there is not way how to store it into data model
    // TODO: extend data model? - or implement using release patterns?
    mcell_warn(
        "Release event %s starts at time different from 0, conversion to data model is not supported yet, ignoring it.",
        release_site_name.c_str()
    );
    return;
  }

  if (release_shape == ReleaseShape::INITIAL_SURF_REGION) {
    // not converting this one - is default
  }
  else if (release_shape == ReleaseShape::LIST) {
    // this release event needs to be split into chunks that use the same
    // species and orientation
    uint release_site_index = 0;
    uint current_index = 0;
    do {
      uint begin_index = current_index;
      species_id_t current_species_id = molecule_list[current_index].species_id;
      orientation_t current_orientation = molecule_list[current_index].orientation;

      // find the largest chunk we can convert as a single release site
      do {
        current_index++;
      } while (
          current_index < molecule_list.size() &&
          molecule_list[current_index].species_id == current_species_id &&
          molecule_list[current_index].orientation == current_orientation
      );

      string name;
      if (begin_index == 0 && current_index == molecule_list.size()) {
        // we are going to generate a single release site
        name = release_site_name;
      }
      else {
        // more release sites - append index
        name = release_site_name + "_" + to_string(release_site_index);
      }

      to_data_model_as_one_release_site(
          mcell_node, current_species_id, current_orientation, name, begin_index, current_index
      );

      release_site_index++;
    } while (current_index < molecule_list.size());
  }
  else {
    // usual case, generate a single release site
    to_data_model_as_one_release_site(
        mcell_node, species_id, orientation, release_site_name, 0, 0
    );
  }
}


bool ReleaseEvent::initialize_walls_for_release() {
  assert(region_expr_root != nullptr);
  assert(release_number_method != ReleaseNumberMethod::Invalid);
  cumm_area_and_pwall_index_pairs.clear();

  // no need to initialize
  const BNG::Species& species = world->get_all_species().get(species_id);
  if (species.is_vol() &&
      (release_number_method == ReleaseNumberMethod::ConstNum ||
       release_number_method == ReleaseNumberMethod::ConcentrationNum
      )
  ) {
    // no need to initialize walls for this case
    return true;
  }

  // only a single region for now
  if (region_expr_root->op != RegionExprOperator::Leaf) {
    return false;
  }

  const Region& reg = world->get_region(region_expr_root->region_id);

  for (auto& it: reg.walls_and_edges) {
    wall_index_t wi = it.first;
    const Wall& w = world->get_partition(PARTITION_ID_INITIAL).get_wall(wi);

    CummAreaPWallIndexPair item;
    item.first = w.area;
    item.second.first = PARTITION_ID_INITIAL;
    item.second.second = wi;

    if (!cumm_area_and_pwall_index_pairs.empty()) {
      item.first += cumm_area_and_pwall_index_pairs.back().first;
    }
    cumm_area_and_pwall_index_pairs.push_back(item);
  }

  return true;
}


static void check_max_release_count(float_t num_to_release, const std::string& name) {
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

    case ReleaseNumberMethod::ConcentrationNum:
      if (diameter == Vec3(LENGTH_INVALID)) {
        // set for instance for ReleaseShape::SPHERICAL
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
            return 0;
            break;

          case ReleaseShape::REGION:
            // number is computed in release_inside_regions
            return 0;

          default:
            mcell_internal_error("Release by concentration on invalid release site "
                                 "shape (%d) for release site \"%s\".",
                                 (int)release_shape, release_site_name.c_str());
            return 0;
            break;
        }
        assert(concentration != FLT_INVALID);
        float_t num_to_release =
            N_AV * 1e-15 * concentration * vol * pow_f(world->config.length_unit, 3.0) + 0.5;
        check_max_release_count(num_to_release, release_site_name);
        return (uint)num_to_release;
      }
      break;

    case ReleaseNumberMethod::DensityNum: {
        // computed in release_onto_regions in MCell3
        assert(!cumm_area_and_pwall_index_pairs.empty());
        float_t max_A = cumm_area_and_pwall_index_pairs.back().first;
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


// returns the number of actually removed molecules
int ReleaseEvent::randomly_remove_molecules(
    Partition& p, const MoleculeIdsVector& mol_ids_in_region, int number_to_remove) {
  // randomly remove molecules
  int num_removed = 0;
  for (size_t i = 0; i < mol_ids_in_region.size(); i++) {
    molecule_id_t m_id = mol_ids_in_region[i];
    int remaining = mol_ids_in_region.size() - i;

    if (rng_dbl(&world->rng) < ((float_t)(number_to_remove)) / ((float_t)remaining)) {
      p.set_molecule_as_defunct(p.get_m(m_id));
      num_removed++;
      number_to_remove--;
    }
  }
  release_assert(number_to_remove >= 0);
  return num_removed;
}


/***************************************************************************
vacuum_from_regions:
  Molecules of the specified type are
       removed uniformly at random from the free area in the regions
       specified by the release site object.
  Note: if the user requests to remove more molecules than actually exist,
        the function will return success and not give a warning.
***************************************************************************/
int ReleaseEvent::vacuum_from_regions(int number_to_remove) {
  assert(!cumm_area_and_pwall_index_pairs.empty());
  assert(number_to_remove > 0);

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  MoleculeIdsVector mol_ids_on_region;

  for (auto& item: cumm_area_and_pwall_index_pairs) {
    assert(item.second.first == PARTITION_ID_INITIAL);
    const Wall& w = p.get_wall(item.second.second);

    for (molecule_id_t m_id: w.grid.get_molecules_per_tile()) {
      if (m_id != MOLECULE_ID_INVALID) {
        const Molecule& m = p.get_m(m_id);
        assert(m.is_surf());
        if (m.is_defunct()) {
          continue;
        }
        if (m.species_id != species_id) {
          continue;
        }
        // NOTE: MCell3 does not care about orientation
        if (orientation != ORIENTATION_NONE && m.s.orientation != orientation) {
          continue;
        }
        mol_ids_on_region.push_back(m_id);
      }
    }
  }

  return randomly_remove_molecules(p, mol_ids_on_region, number_to_remove);
}


void ReleaseEvent::release_onto_regions(int& computed_release_number) {
  int success = 0, failure = 0;
  float_t seek_cost = 0;

  assert(!cumm_area_and_pwall_index_pairs.empty());
  float_t total_area = cumm_area_and_pwall_index_pairs.back().first;
  float_t est_sites_avail = (int)total_area;
  const float_t rel_list_gen_cost = 10.0; /* Just a guess */
  float_t pick_cost = rel_list_gen_cost * est_sites_avail;

  int n = computed_release_number;

  if (n < 0) {
    computed_release_number = -vacuum_from_regions(-n);
    return;
  }

  const int too_many_failures = 10; /* Just a guess */
  while (n > 0) {
    if (failure >= success + too_many_failures) {
      seek_cost =
          n * (((double)(success + failure + 2)) / ((double)(success + 1)));
    }

    if (seek_cost < pick_cost) {
      float_t A = rng_dbl(&world->rng) * total_area;
      size_t cum_area_index = cum_area_bisect_high(cumm_area_and_pwall_index_pairs, A);
      PartitionWallIndexPair pw = cumm_area_and_pwall_index_pairs[cum_area_index].second;
      Partition& p = world->get_partition(pw.first);
      Wall& wall = p.get_wall(pw.second);

      if (!wall.has_initialized_grid()) {
        wall.initialize_grid(p); // sets wall's grid_index
      }

      Grid& grid = wall.grid;

      // get the random number for the current wall
      if (cum_area_index != 0) {
        A -= cumm_area_and_pwall_index_pairs[cum_area_index - 1].first;
      }

      tile_index_t tile_index = (grid.num_tiles_along_axis * grid.num_tiles_along_axis) * (A / wall.area);
      if (tile_index >= grid.num_tiles) {
        tile_index = grid.num_tiles - 1;
      }

      if (grid.get_molecule_on_tile(tile_index) != MOLECULE_ID_INVALID) {
        failure++;
        continue;
      }

      molecule_id_t sm_id =
          GridUtil::place_single_molecule_onto_grid(
              p, world->rng, wall, tile_index, false, Vec2(),
              species_id, orientation, event_time, get_release_delay_time()
          );

      schedule_for_immediate_diffusion_if_needed(sm_id, WallTileIndexPair(wall.index, tile_index));

      #ifdef DEBUG_RELEASES
        p.get_m(sm_id).dump(p, "Released sm:", "", p.stats.get_current_iteration(), actual_release_time, true);
      #endif

      success++;
      n--;
    }
    else {
      // TODO: MCell3 handles these cases better
      const BNG::Species& species = world->get_all_species().get(species_id);
      mcell_error("Could not release %d of %s at %s, too many failed attempts to place surface molecules.",
          computed_release_number, species.name.c_str(), release_site_name.c_str());
    }
  }
}


bool ReleaseEvent::is_point_inside_region_expr_recursively(Partition& p, const Vec3& pos, const RegionExprNode* region_expr_node) {
  assert(region_expr_node->op != RegionExprOperator::Invalid);

  if (region_expr_node->op == RegionExprOperator::Leaf) {
    Region& reg = p.get_region_by_id(region_expr_node->region_id);
    return reg.is_point_inside(p, pos);
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


/*
 * num_vol_mols_from_conc computes the number of volume molecules to be
 * released within a closed object. There are two cases:
 * - for a single closed object we know the exact volume and can thus compute
 *   the exact number of molecules required and release them by setting
 *   exactNumber to true.
 * - for a release object consisting of a boolean expression of closed objects
 *   we are currently not able to compute the volume exactly. Instead we compute
 *   the number of molecules in the bounding box and then release an approximate
 *   number by setting exactNumber to false.
 */
uint ReleaseEvent::num_vol_mols_from_conc(bool &exact_number) {
  Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  float_t vol = 0.0;
  if (region_expr_root->op == RegionExprOperator::Leaf) {
    Region& r = p.get_region_by_id(region_expr_root->region_id);
    r.initialize_volume_info_if_needed(p);
    release_assert(r.is_manifold() && "Trying to release into a regions that is not a manifold and has no volume");
    vol = r.get_volume();
    exact_number = true;
  }
  else {
    // estimate the volume
    vol =
        (region_urb.x - region_llf.x) *
        (region_urb.y - region_llf.y) *
        (region_urb.z - region_llf.z);
    exact_number = false;
  }

  float_t num_to_release =
      N_AV * 1e-15 * concentration * vol * pow_f(world->config.length_unit, 3) + 0.5;

  check_max_release_count(num_to_release, release_site_name);
  return num_to_release;
}


/*************************************************************************
vacuum_inside_regions:
  Note: if more molecules are to be removed than actually exist, all
        existing molecules of the specified type are removed.
  Returns number of actually removed molecules.
*************************************************************************/
int ReleaseEvent::vacuum_inside_regions(int number_to_remove) {
  assert(number_to_remove > 0);

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  MoleculeIdsVector mol_ids_in_region;

  // get all molecules that match the removed species
  // MCell3 knows which molecules are in which subparts, we don't know this
  // so we must go through all molecule
  for (const Molecule& m: p.get_molecules()) {
    if (m.is_defunct()) {
      continue;
    }
    if (m.species_id != species_id) {
      continue;
    }
    release_assert(m.is_vol());

    // filter by bounding box
    if (!point_in_box(m.v.pos, region_llf, region_urb)) {
      continue;
    }
    // then precisely by region
    if (!is_point_inside_region_expr_recursively(p, m.v.pos, region_expr_root)) {
      continue;
    }

    mol_ids_in_region.push_back(m.id);
  }

  return randomly_remove_molecules(p, mol_ids_in_region, number_to_remove);
}


/*************************************************************************
release_inside_regions:
  Note: if the CCNNUM release method is used, the number of molecules
        passed in is ignored.
*************************************************************************/
void ReleaseEvent::release_inside_regions(int& computed_release_number) {

  assert(region_expr_root != nullptr);

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  bool exact_number = false;
  if (release_number_method == ReleaseNumberMethod::ConcentrationNum) {
    computed_release_number = num_vol_mols_from_conc(exact_number);
  }

  int n = computed_release_number;

  if (n < 0) {
    computed_release_number = -vacuum_inside_regions(-n);
    return;
  }

  while (n > 0) {
    Vec3 pos;
    pos.x = region_llf.x + (region_urb.x - region_llf.x) * rng_dbl(&world->rng);
    pos.y = region_llf.y + (region_urb.y - region_llf.y) * rng_dbl(&world->rng);
    pos.z = region_llf.z + (region_urb.z - region_llf.z) * rng_dbl(&world->rng);

    if (!is_point_inside_region_expr_recursively(p, pos, region_expr_root)) {
      if (release_number_method == ReleaseNumberMethod::ConcentrationNum && !exact_number) {
        computed_release_number--;
        n--;
      }
      continue;
    }

    // TODO_LATER: location can be close to a partition boundary, we might need to release to a different partition
    Molecule& new_vm = p.add_volume_molecule(
        Molecule(MOLECULE_ID_INVALID, species_id, pos, event_time), get_release_delay_time()
    );
    new_vm.set_flag(MOLECULE_FLAG_VOL);
    new_vm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

    schedule_for_immediate_diffusion_if_needed(new_vm.id);

    n--;

#ifdef DEBUG_RELEASES
    new_vm.dump(p, "Released vm:", "", p.stats.get_current_iteration(), actual_release_time, true);
#endif
  }
}


void ReleaseEvent::release_ellipsoid_or_rectcuboid(int computed_release_number) {
  assert(computed_release_number >= 0 && "Cannot have negative SPHERICAL release");

  Partition& p = world->get_partition(world->get_or_add_partition_index(location));
  float_t time_step = world->get_all_species().get(species_id).time_step;

  const int is_spheroidal = (release_shape == ReleaseShape::SPHERICAL ||
                             /*release_shape == SHAPE_ELLIPTIC ||*/
                             release_shape == ReleaseShape::SPHERICAL_SHELL);

  for (int i = 0; i < computed_release_number; i++) {
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
        Molecule(MOLECULE_ID_INVALID, species_id, molecule_location, event_time), get_release_delay_time()
    );
    new_vm.set_flag(MOLECULE_FLAG_VOL);
    new_vm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

    schedule_for_immediate_diffusion_if_needed(new_vm.id);

#ifdef DEBUG_RELEASES
    new_vm.dump(p, "Released vm:", "", p.stats.get_current_iteration(), actual_release_time, true);
#endif

  }
}


void ReleaseEvent::release_list() {
  for (const SingleMoleculeReleaseInfo& info: molecule_list) {

    BNG::Species& species = world->get_all_species().get(info.species_id);
    Partition& p = world->get_partition(world->get_or_add_partition_index(info.pos));

    if (species.is_vol()) {
      Molecule& new_vm = p.add_volume_molecule(
          Molecule(MOLECULE_ID_INVALID, info.species_id, info.pos, event_time), get_release_delay_time()
      );
      new_vm.set_flag(MOLECULE_FLAG_VOL);
      new_vm.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

      schedule_for_immediate_diffusion_if_needed(new_vm.id);

      cout
        << "Released 1 " << species.name << " from \"" << release_site_name << "\""
        << " at iteration " << world->get_current_iteration() << ".\n";
    }
    else {
      orientation_t orient;
      assert(info.orientation != ORIENTATION_NOT_SET);
      if (info.orientation == ORIENTATION_NONE) {
        orient = (rng_uint(&world->rng) & 1) ? 1 : -1;
      }
      else {
        orient = info.orientation;
      }

      float_t diam = diameter.x;
      assert(diam != FLT_INVALID);
      molecule_id_t sm_id = GridUtil::place_surface_molecule_to_closest_pos(
          p, world->rng, info.pos, info.species_id, orient, diameter.x,
          event_time, get_release_delay_time()
      );

      const Molecule& sm = p.get_m(sm_id);
      schedule_for_immediate_diffusion_if_needed(sm_id, WallTileIndexPair(sm.s.wall_index, sm.s.grid_tile_index));

      if (sm_id != MOLECULE_ID_INVALID) {
        cout
          << "Released 1 " << species.name << " from \"" << release_site_name << "\""
          << " at iteration " << world->get_current_iteration() << ".\n";
      }
      else {
        // TODO: unify logging, have some C++ streams,,,
        // TODO: copy message from MCell3
        mcell_error("Could not release %s from %s, possibly the release diameter is too short.",
            species.name.c_str(), release_site_name.c_str()
        );
      }

    }
  }
}


void ReleaseEvent::init_surf_mols_by_number(Partition& p, const Region& reg, const InitialRegionMolecules& info) {
  uint n_free_sm = 0;

  /* initialize surface molecule grids in region as needed and */
  /* count total number of free surface molecule sites in region */

  vector<WallTileIndexPair> free_tiles;

  for (auto wall_edge_it: reg.walls_and_edges) {
    Wall& w = p.get_wall(wall_edge_it.first);
    if (!w.has_initialized_grid()) {
      w.initialize_grid(p);
    }

    Grid& g = w.grid;
    n_free_sm += g.get_num_free_tiles();

    for (tile_index_t ti = 0; ti < g.num_tiles; ti++) {
      if (g.get_molecule_on_tile(ti) == MOLECULE_ID_INVALID) {
        free_tiles.push_back(WallTileIndexPair(w.index, ti));
      }
    }
  }
  assert(n_free_sm == free_tiles.size() && "Num free tiles does not match");

  if (info.release_num > 0 && n_free_sm == 0) {
    mcell_error("Number of free surface molecule tiles in region %s = %d", reg.name.c_str(), n_free_sm);
  }

  if (info.release_num > n_free_sm / 2) {
    mcell_warn("Implementation of filling more than half of free tiles is different in MCell4 from MCell3.");
  }

  for (uint i = 0; i < info.release_num; i++) {
    uint num_attempts = 0;
    /* Loop until we find a vacant tile. */
    while (1) {
      uint slot_num = (int)(rng_dbl(&world->rng) * n_free_sm);

      const WallTileIndexPair& wip = free_tiles[slot_num];
      Wall& w = p.get_wall(wip.wall_index);

      if (w.grid.get_molecule_on_tile(wip.tile_index) == MOLECULE_ID_INVALID) {
        GridUtil::place_single_molecule_onto_grid(
            p, world->rng, w, wip.tile_index, false, Vec2(),
            info.species_id, info.orientation, event_time, get_release_delay_time()
        );
        break;
      }

      if (num_attempts != 0 && num_attempts % 10000 == 0) {
        mcell_warn("Made %d of attempts while placing molecule in release event %s.",
            num_attempts, reg.name.c_str());
      }
    }
  }
}


void ReleaseEvent::init_surf_mols_by_density(
    Partition& p, Wall& w,
    map<species_id_t, uint>& num_released_per_species
) {
  if (!w.has_initialized_grid()) {
    w.initialize_grid(p);
  }

  float_t tot_prob = 0;
  float_t tot_density = 0;

  vector<pair<float_t, InitialRegionMolecules>> prob_info_pairs;

  // do for all surface regions of a wall at once
  for (region_index_t reg_index: w.regions) {
    const Region& reg = p.get_region(reg_index);
    for (const InitialRegionMolecules& info: reg.initial_region_molecules) {
      if (!info.is_release_by_density()) {
        // skip release by number, they will be handled later
        continue;
      }

      tot_prob += (w.area * info.release_density) / (w.grid.num_tiles * world->config.grid_density);

      // make an array with cummulative probs
      prob_info_pairs.push_back(make_pair(tot_prob, info));

      tot_density += info.release_density;
    }
  }

  if (tot_density > world->config.grid_density) {
    mcell_warn(
        "Total surface molecule density too high: %f.  Filling all available "
        "surface molecule sites.",
        tot_density);
  }

  if (prob_info_pairs.empty()) {
    // nothing to do
    return;
  }

  // for each tile of the wall
  for (tile_index_t ti = 0; ti < w.grid.num_tiles; ti++) {
    float_t rnd = rng_dbl(&world->rng);

    size_t index;
    for (index = 0; index < prob_info_pairs.size(); ++index) {
      if (rnd <= prob_info_pairs[index].first) {
        break;
      }
    }

    if (index >= prob_info_pairs.size()) {
      continue;
    }

    species_id_t species_id = prob_info_pairs[index].second.species_id;
    GridUtil::place_single_molecule_onto_grid(
        p, world->rng, w, ti, false, Vec2(),
        prob_info_pairs[index].second.species_id, prob_info_pairs[index].second.orientation,
        event_time, get_release_delay_time()
    );

    auto it = num_released_per_species.find(species_id);
    if (it != num_released_per_species.end()) {
      it->second++;
    }
    else {
      num_released_per_species[species_id] = 1;
    }
  }
}


// based on init_wall_surf_mols
void ReleaseEvent::release_initial_molecules_onto_surf_regions() {
  release_assert(running_diffuse_event_to_update == nullptr && "Cannot be executed during diffusion&react event");

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  // first collect all walls that are affected
  vector<wall_index_t> walls;
  for (const Wall& w: p.get_walls()) {
    for (region_index_t reg_index: w.regions) {
      const Region& reg = p.get_region(reg_index);
      if (reg.has_initial_molecules()) {
        walls.push_back(w.index);
        break;
      }
    }
  }

  map<species_id_t, uint> num_released_per_species;
  for (wall_index_t wi: walls) {
    Wall& w = p.get_wall(wi);
    init_surf_mols_by_density(p, w, num_released_per_species);
  }
  for (auto it: num_released_per_species) {
    cout <<
        "Released " << it.second << " " << world->get_all_species().get(it.first).name << " onto surface regions " <<
        "at iteration " << world->get_current_iteration() << " (specified with density).\n";
  }

  for (const Region& reg: p.get_regions()) {
    // for each specifies initial molecules
    for (const InitialRegionMolecules& info: reg.initial_region_molecules) {
      if (!info.is_release_by_num()) {
        // skip density, they were already handled
        continue;
      }
      init_surf_mols_by_number(p, reg, info);

      cout
        << "Released " << info.release_num << " " << world->get_all_species().get(info.species_id).name << " on region \"" << reg.name << "\""
        << " at iteration " << world->get_current_iteration() << " (specified with number).\n";
    }
  }
}

// TODO: cleanup the release number computation
void ReleaseEvent::step() {

  perf() << "Starting release from \"" << release_site_name << "\".\n";

  int num_released = 0;
  if (release_shape == ReleaseShape::REGION) {
    int number = calculate_number_to_release();
    const BNG::Species& species = world->get_all_species().get(species_id);
    if (species.is_surf()) {
      release_onto_regions(number);
    }
    else {
      release_inside_regions(number);
    }
    num_released = number;
  }
  else if (release_shape == ReleaseShape::SPHERICAL) {
    int number = calculate_number_to_release();
    assert(diameter.is_valid());
    release_ellipsoid_or_rectcuboid(number);
    num_released = number;
  }
  else if (release_shape == ReleaseShape::LIST) {
    release_list();
  }
  else if (release_shape == ReleaseShape::INITIAL_SURF_REGION) {
    release_initial_molecules_onto_surf_regions();
  }
  else {
    assert(false);
  }

  if (release_shape == ReleaseShape::REGION || release_shape == ReleaseShape::SPHERICAL) {
    const BNG::Species& species = world->get_all_species().get(species_id);
    const char* type;
    int count;
    if (num_released >= 0) {
      type = "Released ";
      count = num_released;
    }
    else {
      type =  "Removed ";
      count = -num_released;
    }
    cout
      << type << count << " " << species.name << " from \"" << release_site_name << "\""
      << " at iteration " << world->get_current_iteration() << ".\n";
  }
}


void ReleaseEvent::schedule_for_immediate_diffusion_if_needed(
    const molecule_id_t id, const WallTileIndexPair& where_released) {
  // NOTE: we do nto care about partitions here but we should
  if (running_diffuse_event_to_update != nullptr) {
    running_diffuse_event_to_update->add_diffuse_action(DiffuseAction(id, where_released));
  }
}

} // namespace mcell


