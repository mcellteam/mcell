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

#include "config.h"
#include "logging.h"
#include "sym_table.h"
#include "mcell_species.h"
#include "init.h"
#include "mcell_structs.h"

MCELL_STATUS mcell_add_surf_class_properties(
    MCELL_STATE *state,
    int reaction_type,
    mcell_symbol *sc_sym,
    mcell_symbol *reactant_sym,
    short orient);

MCELL_STATUS mcell_create_surf_class(
    MCELL_STATE *state,
    const char *surf_class_name,
    mcell_symbol **sc_sym);

struct sm_dat *mcell_add_mol_release_to_surf_class(
    MCELL_STATE *state,
    struct sym_entry *sc_sym,
    struct mcell_species *sm_info,
    double quantity,
    int density_or_num,
    struct sm_dat *smd_list);

MCELL_STATUS mcell_assign_surf_class_to_region(
    struct sym_entry *sc_sym,
    struct region *rgn);
