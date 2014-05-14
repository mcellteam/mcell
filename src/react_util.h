/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2014 by *
 * The Salk Institute for Biological Studies and *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or *
 * modify it under the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 *
 * of the License, or (at your option) any later version. *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *
 * GNU General Public License for more details. *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License *
 * along with this program; if not, write to the Free Software *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 *USA. *
 *                                                                                 *
 ***********************************************************************************/

#ifndef REACT_UTIL_H
#define REACT_UTIL_H

#include "mcell_structs.h"

double compute_pb_factor(struct volume *world, struct rxn *rx,
                         int max_num_surf_products);

int get_rxn_by_name(struct rxn **reaction_hash, int hashsize,
                    const char *rx_name, struct rxn **found_rx, int *path_id);

int change_reaction_probability(struct volume *world, struct rxn *rx,
                                int path_id, double new_rate);

void issue_reaction_probability_warnings(struct volume *world, struct rxn *rx);

#endif
