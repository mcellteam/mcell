/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by *
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

#ifndef CREATE_SPECIES_H
#define CREATE_SPECIES_H
#include "libmcell.h"

struct species *assemble_mol_species(MCELL_STATE *state,
                                     struct sym_table *sym_ptr,
                                     struct mcell_species_spec *species);

int new_mol_species(MCELL_STATE *state, char *name, struct sym_table *sym_ptr);

int ensure_rdstep_tables_built(MCELL_STATE *state);

#endif
