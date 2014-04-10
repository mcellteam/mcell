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
  struct mcell_species *reactants, int *num_reactants, int *num_vol_mols,
  int *num_grid_mols, int *num_complex_reactants, int *all_3d, 
  int *oriented_count, int *complex_type);

MCELL_STATUS extract_catalytic_arrow(struct pathway *path, 
  struct reaction_arrow *react_arrow, int *num_reactants, 
  int *num_vol_mols, int *num_grid_mols, int *all_3d,
  int *oriented_count);

MCELL_STATUS extract_surface(struct pathway *path, 
  struct mcell_species *surf_class, int *num_reactants, 
  unsigned int *num_surfaces, int *oriented_count);

MCELL_STATUS extract_products(MCELL_STATE *state, struct pathway *path, 
  struct mcell_species *products, int *num_surf_products, 
  int *num_complex_products, int bidirectional, int complex_type,
  int all_3d);

MCELL_STATUS check_surface_specs(MCELL_STATE *state, int num_reactants, 
    int num_surfaces, int num_vol_mols, int all_3d, int oriented_count);

MCELL_STATUS extract_pathname(struct pathway *path, struct rxn *rxnp, 
    struct sym_table *pathname);

MCELL_STATUS add_catalytic_species_to_products(struct pathway *path, 
    int catalytic, int bidirectional, int all_3d);

MCELL_STATUS create_product_signature(struct pathway *path);

MCELL_STATUS extract_forward_rate(struct pathway *path, 
    struct reaction_rates *rate, const char *rate_filename);

MCELL_STATUS grid_space_available_for_surface_products(
    double vacancy_search_dist2, int num_grid_mols, int num_vol_mols, 
    int num_surf_products);

MCELL_STATUS invert_current_reaction_pathway(MCELL_STATE *state,
    struct pathway *pathp, struct reaction_rate *reverse_rate,
    const char *rate_filename);

char* create_rx_name(struct pathway *p);
char* create_prod_signature(struct product **product_head);

MCELL_STATUS init_reactions(MCELL_STATE *state);

#endif
