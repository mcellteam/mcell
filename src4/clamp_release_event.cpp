/******************************************************************************
 *
 * Copyright (C) 2020 by
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

#include "clamp_release_event.h"

#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs.h"
#include "logging.h"

#include <iostream>

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
  cout << ind2 << "concentration: \t\t" << concentration << " [float_t]\n";
  cout << ind2 << "scaling_factor: \t\t" << scaling_factor << " [float_t]\n";
  dump_cumm_area_and_pwall_index_pairs(cumm_area_and_pwall_index_pairs, ind2);
}


void ClampReleaseEvent::to_data_model(Json::Value& mcell_node) const {
  // --- define_surface_classes --- (only when needed)
  Json::Value& define_surface_classes = mcell_node[KEY_DEFINE_SURFACE_CLASSES];
  DMUtil::add_version(define_surface_classes, VER_DM_2014_10_24_1638);
  Json::Value& surface_class_list = define_surface_classes[KEY_SURFACE_CLASS_LIST];

  Json::Value clamp;
  clamp[KEY_NAME] = world->get_all_species().get(surf_class_species_id).name;
  clamp[KEY_DESCRIPTION] = "";
  DMUtil::add_version(clamp, VER_DM_2018_01_11_1330);

  Json::Value& surface_prop_list = clamp[KEY_SURFACE_CLASS_PROP_LIST];
  surface_prop_list = Json::Value(Json::arrayValue);

  Json::Value surface_prop_item;
  surface_prop_item[KEY_CLAMP_VALUE] = f_to_str(concentration);
  surface_prop_item[KEY_SURF_CLASS_ORIENT] = DMUtil::orientation_to_str(orientation);
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
  DMUtil::add_version(surface_prop_item, VER_DM_2015_11_08_1756);

  surface_prop_list.append(surface_prop_item);
  surface_class_list.append(clamp);

  // the affected region is printed in Partition::to_data_model
}


void ClampReleaseEvent::update_cumm_areas_and_scaling() {
  float_t cumm_area = 0;
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
static int poisson_dist(float_t lambda, float_t p) {
  int i, lo, hi;
  float_t plo, phi, pctr;
  float_t lambda_i;

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

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  const BNG::Species& species = world->bng_engine.get_all_species().get(species_id);
  assert(species.is_vol());
  float_t n_collisions = scaling_factor * species.space_step * concentration / species.time_step;
  if (orientation != ORIENTATION_NONE) {
    n_collisions *= 0.5;
  }
  int n_to_emit = poisson_dist(n_collisions, rng_dbl(&world->rng));
  if (n_to_emit == 0) {
    return;
  }

  float_t total_area = cumm_area_and_pwall_index_pairs.back().first;

  while (n_to_emit > 0) {

    int idx = cum_area_bisect_high(cumm_area_and_pwall_index_pairs, rng_dbl(&world->rng) * total_area);

    assert(cumm_area_and_pwall_index_pairs[idx].second.first == PARTITION_ID_INITIAL);
    const Wall& w = p.get_wall(cumm_area_and_pwall_index_pairs[idx].second.second);
    const Vec3& v0 = p.get_wall_vertex(w, 0);
    const Vec3& v1 = p.get_wall_vertex(w, 1);
    const Vec3& v2 = p.get_wall_vertex(w, 2);

    float_t s1 = sqrt_f(rng_dbl(&world->rng));
    float_t s2 = rng_dbl(&world->rng) * s1;

    Vec3 v;
    v = v0 + Vec3(s1) * (v1 - v0) + Vec3(s2) * (v2 - v1);

    orientation_t o;
    if (orientation == ORIENTATION_NONE) {
      o = (rng_uint(&world->rng) & 2) - 1;
    }
    else {
      o = orientation;
    }

    float_t eps = EPS * o;

    float_t move1 = abs_f(v.x);
    float_t move2 = abs_f(v.y);
    if (move1 < move2) {
      move1 = move2;
    }
    move2 = abs_f(v.z);
    if (move1 < move2) {
      move1 = move2;
    }
    if (move1 > 1.0){
      eps *= move1;
    }

    Vec3 pos = v + w.normal * Vec3(eps);

    Molecule& new_vm = p.add_volume_molecule(
          Molecule(MOLECULE_ID_INVALID, species_id, pos, event_time), 0.5 /* release delay time */
    );
    new_vm.flags |= ACT_CLAMPED | IN_VOLUME | ACT_DIFFUSE | MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN;

    new_vm.set_clamp_orientation(o);
    new_vm.v.previous_wall_index = w.index;

    n_to_emit--;
  }
}

} // namespace mcell


