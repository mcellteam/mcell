/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
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

#ifndef REACT_NFSIM_H
#define REACT_NFSIM_H

#include "mcell_structs.h"
 
/*
calculates particle orientation based on nfsim compartment information
*/
void calculate_reactant_orientation(struct abstract_molecule* reac, struct abstract_molecule* reac2, 
                            bool* orientation_flag1, bool* orientation_flag2, 
                            int* reactantOrientation1, int* reactantOrientation2);

queryOptions initializeNFSimQueryForBimolecularReactions(struct graph_data *am, 
                                                      struct graph_data* am2,
                                                      char* onlyActive);

int trigger_bimolecular_preliminary_nfsim(struct abstract_molecule *reacA,
                                    struct abstract_molecule *reacB);

void pick_unimolecular_reaction_nfsim(struct volume *state,
                                       struct abstract_molecule *am, struct rxn* rx);

#endif