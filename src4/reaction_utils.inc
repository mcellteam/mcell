/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 *
 * There is an exception, methods that schedule actions into the event's
 * new_diffuse_or_unimol_react_actions action queue stay in the diffuse_react_event_t
 * class.
 */

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

using namespace std;

namespace mcell {

namespace rx_util {

// ---------------------------------- bimolecular reactions ----------------------------------

/*************************************************************************
test_bimolecular
  In: the reaction we're testing
      a scaling coefficient depending on how many timesteps we've
        moved at once (1.0 means one timestep) and/or missing interaction area
      local probability factor (positive only for the reaction between two
        surface molecules, otherwise equal to zero)
      reaction partners
  Out: RX_NO_RX if no reaction occurs
       int containing which reaction pathway to take if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
*************************************************************************/
/*, double scaling - 1, double local_prob_factor = 0,*/
static int test_bimolecular(
    const reaction_t& rx,
    rng_state& rng,
    const volume_molecule_t& a1, // unused for now
    const volume_molecule_t& a2, // unused for now
    const float_t scaling
) {
  /* rescale probabilities for the case of the reaction
     between two surface molecules */
  float_t min_noreaction_p = rx.min_noreaction_p; // local_prob_factor == 0

  assert(min_noreaction_p < scaling);
  /* Instead of scaling rx->cum_probs array we scale random probability */
  float_t p = rng_dbl(&rng) * scaling;

  if (p >= min_noreaction_p) {
    return RX_NO_RX;
  }
  else {
    return 0; // we have just one pathwayy
  }
}


// might return nullptr if there is no unimolecular reaction for this species
// based on pick_unimolecular_reaction
static const reaction_t* pick_unimol_rx(
    const world_t* world,
    const species_id_t species_id
) {
  const unimolecular_reactions_map_t* unimol_rxs = world->world_constants.unimolecular_reactions_map;
  auto it = unimol_rxs->find(species_id);
  if (it == unimol_rxs->end()) {
    return nullptr;
  }
  else {
    return it->second;
  }
}


// based on timeof_unimolecular
static float_t time_of_unimol(const reaction_t* rx, rng_state& rng) {
  double k_tot = rx->max_fixed_p;
  double p = rng_dbl(&rng);

  if ((k_tot <= 0) || (!distinguishable(p, 0, EPS)))
    return TIME_FOREVER;
  return -log(p) / k_tot;
}


// based on compute_lifetime
static float_t compute_unimol_lifetime(
    world_t* world,
    const volume_molecule_t& vm,
    const reaction_t* rx
) {
  assert(rx != nullptr);

  float_t res = time_of_unimol(rx, world->rng);

#ifdef DEBUG_REACTIONS
  DUMP_CONDITION4(
      // calling rng for unimolecular
      vm.dump(world, "Assigned unimolecular time (prev rng):", "", world->current_iteration, res);
  );
#endif

  return res;
}

} // namespace rx_util

} // namespace mcell
