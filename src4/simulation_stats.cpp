/******************************************************************************
 *
 * Copyright (C) 2019-2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "simulation_stats.h"

#include "bng/bng.h"

using namespace std;

namespace MCell {


RxnCountStats& SimulationStats::get_rxn_stats(
    const BNG::RxnContainer& all_rxns, const BNG::RxnClass* rxn_class) {
  // we must use reactant name pairs to store the information
  // because their count does not grow and rxn classes

  // also assuming that the names of reactants in reactions are always sorted
  // (does not matter how), i.e. there is always A + B -> C and A + B -> D,
  // and not B + A -> E

  assert(rxn_class->get_num_reactions() >= 1);
  BNG::rxn_rule_id_t id = rxn_class->get_rxn_rule_id(0);

  const BNG::RxnRule* rxn = all_rxns.get(id);
  assert(rxn != nullptr);
  assert(rxn->is_bimol());

  const string& n1 = rxn->reactants[0].name;
  assert(n1 != "");
  const string& n2 = rxn->reactants[1].name;
  assert(n2 != "");

  // increment or add a new item if the pair we are searching for does not exist
  auto it1 = bimol_rxn_stats.find(n1);
  if (it1 == bimol_rxn_stats.end()) {
    it1 = bimol_rxn_stats.insert(make_pair(n1, NameStatsMap())).first;
  }

  auto it2 = it1->second.find(n2);
  if (it2 == it1->second.end()) {
    it2 = it1->second.insert(make_pair(n2, RxnCountStats())).first;
  }

  return it2->second;
}

void SimulationStats::inc_rxn_occurred(
    const BNG::RxnContainer& all_rxns, const BNG::RxnClass* rxn_class,
    const uint64_t occurred) {

  get_rxn_stats(all_rxns, rxn_class).occurred += occurred;
}


void SimulationStats::inc_rxn_skipped(
    const BNG::RxnContainer& all_rxns, const BNG::RxnClass* rxn_class,
    const double skipped) {

  get_rxn_stats(all_rxns, rxn_class).skipped += skipped;
}

void SimulationStats::print_missed_rxns_warnings() {
  bool first_report = true;
  for (auto n1_map_pair: bimol_rxn_stats) {
    for (auto n2_stats_pair: n1_map_pair.second) {
      const RxnCountStats& stats = n2_stats_pair.second;
      if (stats.skipped > 0) {
        if (first_report) {
          cerr << "Warning: Some reactions were missed because total reaction probability exceeded 1.\n";
          first_report = false;
        }

        double perc =
            0.001 * round(1000 * stats.skipped * 100 / (stats.skipped + stats.occurred));

        if (perc >= SQRT_EPS) { // ignore really small values
          cerr <<
              "  " << n1_map_pair.first << " + " << n2_stats_pair.first << "  --  " <<
              perc << "% of reactions missed.\n";
        }
      }
    }
  }
}


void SimulationStats::print_report() {

  cout << "Total number of ray-subvolume intersection tests (number of ray_trace calls): " << ray_voxel_tests << "\n";
  cout << "Total number of ray-polygon intersection tests: " << ray_polygon_tests << "\n";
  cout << "Total number of ray-polygon intersections: " << ray_polygon_colls << "\n";
  cout << "Total number of mol reflections from a wall: " << mol_wall_reflections << "\n";
  cout << "Total number of vol mol vol mol collisions: " << vol_mol_vol_mol_collisions << "\n";
  cout << "Total number of molecule moves between walls: " << mol_moves_between_walls << "\n";
  cout << "Total number of usages of waypoints for counted volumes: " << num_waypoints_used << "\n";
  cout << "Total number of counted volume recomputations: " << recomputations_of_counted_volume << "\n";
  cout << "Total number of diffuse 3d calls: " << diffuse_3d_calls << "\n";
  if (diffusion_number != 0) {
    cout << "Average diffusion jump was: " <<
        diffusion_cummtime / (double)diffusion_number << " timesteps " <<
        " (" << diffusion_cummtime << "/" << diffusion_number << ")\n";
  }

  print_missed_rxns_warnings();
}


} // namespace MCell
