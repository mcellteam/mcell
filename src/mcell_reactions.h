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

#ifndef MCELL_REACTIONS_H
#define MCELL_REACTIONS_H

#include "mcell_species.h"

#define REGULAR_ARROW 0x00
#define ARROW_BIDIRECTIONAL 0x01
#define ARROW_CATALYTIC 0x02

// typedef struct sym_table mcell_symbol;

enum {
  RATE_UNSET = -1,
  RATE_CONSTANT = 0,
  RATE_FILE = 1,
  RATE_COMPLEX = 2
};

/* Special pathway types. */
enum special_pathway_t {
  RFLCT,  /* Special pathway: reflective surface */
  TRANSP, /* Special pathway: transparent surface */
  SINK    /* Special pathway: absorptive surface */
};

struct reaction_def {
  struct sym_table *sym;
};

struct release_single_molecule_list {
  struct release_single_molecule *rsm_head;
  struct release_single_molecule *rsm_tail;
  int rsm_count;
};

struct reaction_arrow {
  int flags;
  struct mcell_species catalyst;
};

struct reaction_rate {
  int rate_type;
  union {
    double rate_constant;
    char *rate_file;
    struct complex_rate *rate_complex;
  } v;
};

struct reaction_rates {
  struct reaction_rate forward_rate;
  struct reaction_rate backward_rate;
};

MCELL_STATUS
mcell_add_reaction(struct notifications *notify,
                   double **r_step_release,
                   struct sym_table_head *rxn_sym_table,
                   u_int radial_subdivisions,
                   double vacancy_search_dist2,
                   struct mcell_species *reactants,
                   struct reaction_arrow *react_arrow,
                   struct mcell_species *surf_class,
                   struct mcell_species *products, struct sym_table *pathname,
                   struct reaction_rates *rates, const char *rate_filename);

MCELL_STATUS mcell_add_surface_reaction(struct sym_table_head *rxn_sym_table,
                                        int reaction_type,
                                        struct species *surface_class,
                                        struct sym_table *reactant_sym,
                                        short orient);

MCELL_STATUS
mcell_add_concentration_clamp(struct sym_table_head *rxn_sym_table,
                              struct species *surface_class,
                              struct sym_table *mol_sym, short orient,
                              double conc);

MCELL_STATUS init_reactions(MCELL_STATE *state);

MCELL_STATUS mcell_change_reaction_rate(MCELL_STATE *state,
                                        const char *reaction_name,
                                        double new_rate);

struct reaction_rates mcell_create_reaction_rates(int forwardRateType,
                                                  int forwardRate,
                                                  int backwardRateType,
                                                  int backwardRate);

//helper functions for initializing a reaction structure

MCELL_STATUS extract_reactants(struct pathway *path,
                                      struct mcell_species *reactants,
                                      int *num_reactants, int *num_vol_mols,
                                      int *num_surface_mols,
                                      int *num_complex_reactants, int *all_3d,
                                      int *oriented_count, int *complex_type);


MCELL_STATUS extract_products(struct notifications *notify,
                                     struct pathway *path,
                                     struct mcell_species *products,
                                     int *num_surf_products,
                                     int *num_complex_products,
                                     int bidirectional, int complex_type,
                                     int all_3d);

int set_product_geometries(struct pathway *path, struct rxn *rx,
                                  struct product *prod);

//initializes the rxn->info field with null values
int init_reaction_info(struct rxn* rx);

char *create_rx_name(struct pathway *p);
#endif
