/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef REACT_NFSIM_H
#define REACT_NFSIM_H

#include "mcell_structs.h"

/*
calculates particle orientation based on nfsim compartment information
*/
void calculate_reactant_orientation(struct abstract_molecule *reac,
                                    struct abstract_molecule *reac2,
                                    bool *orientation_flag1,
                                    bool *orientation_flag2,
                                    int *reactantOrientation1,
                                    int *reactantOrientation2);

queryOptions initializeNFSimQueryForBimolecularReactions(struct graph_data *am,
                                                         struct graph_data *am2,
                                                         const char *onlyActive);

int trigger_bimolecular_preliminary_nfsim(struct abstract_molecule *reacA,
                                          struct abstract_molecule *reacB);

struct rxn *pick_unimolecular_reaction_nfsim(struct volume *state,
                                             struct abstract_molecule *am);

#endif
