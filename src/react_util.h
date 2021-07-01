/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#pragma once

#include "mcell_structs.h"

double compute_pb_factor(double time_unit,
                         double length_unit,
                         double grid_density,
                         double rx_radius_3d,
                         struct reaction_flags *rxn_flags,
                         int *create_shared_walls_info_flag,
                         struct rxn *rx,
                         int max_num_surf_products);

int get_rxn_by_name(struct rxn **reaction_hash, int hashsize,
                    const char *rx_name, struct rxn **found_rx, int *path_id);

int change_reaction_probability(byte *reaction_prob_limit_flag,
                                struct notifications *notify, struct rxn *rx,
                                int path_id, double new_rate);

void issue_reaction_probability_warnings(struct notifications *notify,
                                         struct rxn *rx);
