/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

#ifndef CREATE_REACTIONS_H
#define CREATE_REACTIONS_H
#include "libmcell.h"


MCELL_STATUS extract_reactants(struct pathway *path, 
  struct species_opt_orient *reactants, int *num_reactants, int *num_vol_mols,
  int *num_grid_mols);

MCELL_STATUS extract_catalytic_arrow(struct pathway *path, 
  struct reaction_arrow *react_arrow, int *num_reactants, 
  int *num_vol_mols, int *num_grid_mols);

MCELL_STATUS extract_surface(struct pathway *path, 
    struct species_opt_orient *surf_class, int *num_reactants, 
    int *num_surfaces, int *oriented_count);

#endif
