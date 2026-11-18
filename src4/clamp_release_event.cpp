/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "clamp_release_event.h"

#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs_shared.h"
#include "logging.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include "defines.h"
#include "world.h"
#include "partition.h"
#include "datamodel_defines.h"
#include "release_event.h"


using namespace std;

namespace MCell {

void ClampReleaseEvent::dump(const std::string indent) const {
  cout << "Concentration clamp event:\n";
  string ind2 = indent + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "species_id: \t\t" << species_id << " [species_id_t]\n";
  cout << ind2 << "surf_class_species_id: \t\t" << surf_class_species_id << " [species_id_t]\n";
  cout << ind2 << "concentration: \t\t" << concentration << " [double]\n";
  cout << ind2 << "scaling_factor: \t\t" << scaling_factor << " [double]\n";
  dump_cumm_area_and_pwall_index_pairs(cumm_area_and_pwall_index_pairs, ind2);
}


void ClampReleaseEvent::to_data_model(Json::Value& mcell_node) const {
  // --- define_surface_classes --- (only when needed)
  Json::Value& define_surface_classes = mcell_node[KEY_DEFINE_SURFACE_CLASSES];
  DMUtils::add_version(define_surface_classes, VER_DM_2014_10_24_1638);
  Json::Value& surface_class_list = define_surface_classes[KEY_SURFACE_CLASS_LIST];

  Json::Value clamp;
  clamp[KEY_NAME] = world->get_all_species().get(surf_class_species_id).name;
  clamp[KEY_DESCRIPTION] = "";
  DMUtils::add_version(clamp, VER_DM_2018_01_11_1330);

  Json::Value& surface_prop_list = clamp[KEY_SURFACE_CLASS_PROP_LIST];
  surface_prop_list = Json::Value(Json::arrayValue);

  Json::Value surface_prop_item;
  surface_prop_item[KEY_CLAMP_VALUE] = f_to_str(concentration);
  surface_prop_item[KEY_SURF_CLASS_ORIENT] = DMUtils::orientation_to_str(orientation);
  surface_prop_item[KEY_MOLECULE] = world->get_all_species().get(species_id).name;
  surface_prop_item[KEY_NAME] = ""; // blender exports name but it does not seem to be needed
  surface_prop_item[KEY_AFFECTED_MOLS] = VALUE_SINGLE;
  if (type == ClampType::CONCENTRATION) {
    surface_prop_item[KEY_SURF_CLASS_TYPE] = VALUE_CLAMP_CONCENTRATION;
  }
  else if (type == ClampType::FLUX) {
    surface_prop_item[KEY_SURF_CLASS_TYPE] = VALUE_CLAMP_FLUX;
  }
  else {
    assert(false);
  }
  DMUtils::add_version(surface_prop_item, VER_DM_2015_11_08_1756);

  surface_prop_list.append(surface_prop_item);
  surface_class_list.append(clamp);

  // the affected region is printed in Partition::to_data_model
}


void ClampReleaseEvent::update_cumm_areas_and_scaling() {
  pos_t cumm_area = 0;
  for (auto& area_pindex_pair: cumm_area_and_pwall_index_pairs) {
    PartitionWallIndexPair pindex = area_pindex_pair.second;
    const Wall& w = world->get_partition(pindex.first).get_wall(pindex.second);
    assert(w.area != 0);
    cumm_area += w.area;
    area_pindex_pair.first = cumm_area;
  }

  if (!cumm_area_and_pwall_index_pairs.empty()) {
    scaling_factor =
        cumm_area_and_pwall_index_pairs.back().first * pow_f(world->config.length_unit, 3)
        / 2.9432976599069717358e-9; /* sqrt(MY_PI)/(1e-15*N_AV) */
  }
  else {
    scaling_factor = 0;
  }
}



/*************************************************************************
poisson_dist:
  In: mean value
      random number distributed uniformly between 0 and 1
  Out: integer sampled from the Poisson distribution.
  Note: This does not sample the CDF.  Instead, it works its way outwards
        from the peak of the PDF.  Kinda weird.  It is not the case
        that low values of the random number will give low values.  It
        is also not super-efficient, but it works.
*************************************************************************/
static int poisson_dist(double lambda, double p) {
  int i, lo, hi;
  double plo, phi, pctr;
  double lambda_i;

  i = (int)lambda;
  pctr = exp_f(-lambda + i * log_f(lambda) -
             lgamma(i + 1)); /* Highest probability bin */

  if (p < pctr) {
    return i;
  }

  lo = hi = i;
  plo = phi = pctr; /* Start at highest-probability bin and work outwards */

  p -= pctr;
  lambda_i = 1.0 / lambda;
  while (p > 0) { /* Keep going until we exhaust probabilities */

    if (lo > 0) { /* We still have a low tail, test it */
      plo *= lo * lambda_i; /* Recursive formula for p for this bin */
      lo--;
      if (p < plo)
        return lo;
      p -= plo;
    }
    /* Always test the high tail (it's infinite) */
    hi++;
    phi = phi * lambda / hi; /* Recursive formula for p for this bin */
    if (p < phi) {
      return hi;
    }
    p -= phi + DBL_EPSILON; /* Avoid infinite loop from poor roundoff */
  }

  /* should never get here */
  assert(false);
  return -1;
}


/*************************************************************************
run_concentration_clamp:
  Out: No return value.  Molecules are released at concentration-clamped
       surfaces to maintain the desired concentation.
*************************************************************************/
void ClampReleaseEvent::step() {

  if (cumm_area_and_pwall_index_pairs.empty()) {
    // nothing to do - maybe not schedule at all?
    return;
  }

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  const BNG::Species& species = world->bng_engine.get_all_species().get(species_id);
  assert(species.is_vol());
  double n_collisions = scaling_factor * species.space_step * concentration / species.time_step;
  release_assert(n_collisions < (double)INT_MAX);
  if (orientation != ORIENTATION_NONE) {
    n_collisions *= 0.5;
  }
  int n_to_emit = poisson_dist(n_collisions, rng_dbl(&world->rng));
  if (n_to_emit == 0) {
    return;
  }

  pos_t total_area = cumm_area_and_pwall_index_pairs.back().first;

  while (n_to_emit > 0) {

    int idx = cum_area_bisect_high(cumm_area_and_pwall_index_pairs, rng_dbl(&world->rng) * total_area);

    assert(cumm_area_and_pwall_index_pairs[idx].second.first == PARTITION_ID_INITIAL);
    const Wall& w = p.get_wall(cumm_area_and_pwall_index_pairs[idx].second.second);
    const Vec3& v0 = p.get_wall_vertex(w, 0);
    const Vec3& v1 = p.get_wall_vertex(w, 1);
    const Vec3& v2 = p.get_wall_vertex(w, 2);

    double s1 = sqrt_f(rng_dbl(&world->rng));
    double s2 = rng_dbl(&world->rng) * s1;

    Vec3 v;
    v = v0 + Vec3(s1) * (v1 - v0) + Vec3(s2) * (v2 - v1);

    orientation_t o;
    if (orientation == ORIENTATION_NONE) {
      o = (rng_uint(&world->rng) & 2) - 1;
    }
    else {
      o = orientation;
    }

    pos_t eps = POS_EPS * o;

    pos_t move1 = fabs_p(v.x);
    pos_t move2 = fabs_p(v.y);
    if (move1 < move2) {
      move1 = move2;
    }
    move2 = fabs_p(v.z);
    if (move1 < move2) {
      move1 = move2;
    }
    if (move1 > 1.0f){
      eps *= move1;
    }

    Vec3 pos = v + w.normal * Vec3(eps);

    Molecule& new_vm = p.add_volume_molecule(
          Molecule(MOLECULE_ID_INVALID, species_id, pos, event_time), 0.5 /* release delay time */
    );
    new_vm.flags |= MOLECULE_FLAG_ACT_CLAMPED | MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN;

    new_vm.set_clamp_orientation(o);
    new_vm.v.previous_wall_index = w.index;

    n_to_emit--;
  }
}


// the part of a region's name after the object's name, e.g. "top" for
// "Scene.box1,top"
static std::string region_name_in_object(const std::string& name) {
  const size_t comma = name.rfind(',');
  return (comma == std::string::npos) ? name : name.substr(comma + 1);
}

// "a", "a and b", "a, b and c"
static std::string join_with_and(const std::vector<std::string>& items) {
  std::string res;
  for (size_t i = 0; i < items.size(); i++) {
    if (i != 0) {
      res += (i + 1 == items.size()) ? " and " : ", ";
    }
    res += items[i];
  }
  return res;
}

// a clamp on one side of a wall and a clamp on another act on the same side
// when the sides are equal or either clamps both sides
static bool clamp_sides_meet(const orientation_t a, const orientation_t b) {
  return a == b || a == ORIENTATION_NONE || b == ORIENTATION_NONE;
}


std::vector<std::string> create_clamp_release_events(World* world, const std::vector<ClampSpec>& clamps_as_given) {
  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const std::vector<Region>& regions = p.get_regions();
  const BNG::SpeciesContainer& all_species = world->get_all_species();

  // one order for MDL files and Python API models alike, whatever order the
  // converter gives: by surface class name, then molecule name, then side. It
  // is the order of the events, and so of their random numbers
  std::vector<ClampSpec> clamps = clamps_as_given;
  std::stable_sort(clamps.begin(), clamps.end(),
      [&all_species](const ClampSpec& a, const ClampSpec& b) {
        const std::string& class_a = all_species.get(a.surf_class_species_id).name;
        const std::string& class_b = all_species.get(b.surf_class_species_id).name;
        if (class_a != class_b) {
          return class_a < class_b;
        }
        const std::string& molecule_a = all_species.get(a.species_id).name;
        const std::string& molecule_b = all_species.get(b.species_id).name;
        if (molecule_a != molecule_b) {
          return molecule_a < molecule_b;
        }
        return a.orientation < b.orientation;
      });

  // each clamp brought to a wall by a region of the clamp's surface class,
  // for every wall, in wall order
  struct Application {
    size_t clamp;
    region_index_t region;
  };
  std::map<wall_index_t, std::vector<Application>> applications_of_wall;
  for (size_t c = 0; c < clamps.size(); c++) {
    for (region_index_t ri = 0; ri < regions.size(); ri++) {
      if (regions[ri].species_id != clamps[c].surf_class_species_id) {
        continue;
      }
      for (const auto& wall_and_edges: regions[ri].walls_and_edges) {
        applications_of_wall[wall_and_edges.first].push_back(Application{c, ri});
      }
    }
  }

  // walls clamped more than once for the same volume molecule on the same
  // side, by object and molecule: the walls and, per clamp, the regions that
  // bring it to them
  struct Twice {
    std::set<wall_index_t> walls;
    std::map<size_t, std::set<region_index_t>> regions_of_clamp;
  };
  std::map<std::pair<geometry_object_index_t, species_id_t>, Twice> clamped_twice;
  for (const auto& wall_and_applications: applications_of_wall) {
    const wall_index_t wi = wall_and_applications.first;
    const std::vector<Application>& apps = wall_and_applications.second;
    for (size_t i = 0; i < apps.size(); i++) {
      for (size_t j = i + 1; j < apps.size(); j++) {
        const ClampSpec& a = clamps[apps[i].clamp];
        const ClampSpec& b = clamps[apps[j].clamp];
        if (a.species_id != b.species_id || !clamp_sides_meet(a.orientation, b.orientation)) {
          continue;
        }
        Twice& twice = clamped_twice[std::make_pair(p.get_wall(wi).object_index, a.species_id)];
        twice.walls.insert(wi);
        twice.regions_of_clamp[apps[i].clamp].insert(apps[i].region);
        twice.regions_of_clamp[apps[j].clamp].insert(apps[j].region);
      }
    }
  }

  if (!clamped_twice.empty()) {
    std::vector<std::string> lines;
    for (const auto& key_and_twice: clamped_twice) {
      const Twice& twice = key_and_twice.second;
      // the clamps in their order, the regions by name
      std::vector<std::string> by_class;
      for (const auto& clamp_and_regions: twice.regions_of_clamp) {
        std::set<std::string> names;
        for (const region_index_t ri: clamp_and_regions.second) {
          names.insert(region_name_in_object(regions[ri].name));
        }
        const std::vector<std::string> region_names(names.begin(), names.end());
        by_class.push_back(
            all_species.get(clamps[clamp_and_regions.first].surf_class_species_id).name +
            " through the region" + (region_names.size() == 1 ? " " : "s ") + join_with_and(region_names));
      }
      const size_t n = twice.walls.size();
      lines.push_back(
          "[CLAMP_ON_WALL_TWICE] " + p.get_geometry_object(key_and_twice.first.first).name + ": " +
          std::to_string(n) + (n == 1 ? " wall is" : " walls are") +
          " clamped more than once for the volume molecule " +
          all_species.get(key_and_twice.first.second).name + " (surface class" +
          (by_class.size() == 1 ? " " : "es ") + join_with_and(by_class) + ")");
    }
    return lines;
  }

  // one event per clamp, over every wall that carries its surface class
  for (size_t c = 0; c < clamps.size(); c++) {
    ClampReleaseEvent* clamp_event = new ClampReleaseEvent(world);
    clamp_event->type = clamps[c].type;
    clamp_event->species_id = clamps[c].species_id;
    clamp_event->surf_class_species_id = clamps[c].surf_class_species_id;
    clamp_event->orientation = clamps[c].orientation;
    clamp_event->concentration = clamps[c].concentration;

    // run each timestep
    clamp_event->event_time = 0;
    clamp_event->periodicity_interval = 1;

    for (const auto& wall_and_applications: applications_of_wall) {
      for (const Application& app: wall_and_applications.second) {
        if (app.clamp == c) {
          clamp_event->cumm_area_and_pwall_index_pairs.push_back(
              CummAreaPWallIndexPair(0, PartitionWallIndexPair(PARTITION_ID_INITIAL, wall_and_applications.first)));
          break;
        }
      }
    }
    clamp_event->update_cumm_areas_and_scaling();
    world->scheduler.schedule_event(clamp_event);
  }
  return std::vector<std::string>();
}

} // namespace mcell
