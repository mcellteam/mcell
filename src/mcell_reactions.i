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

#define REGULAR_ARROW 0x00
#define ARROW_BIDIRECTIONAL 0x01
#define ARROW_CATALYTIC 0x02

// typedef struct sym_entry mcell_symbol;

enum {
  RATE_UNSET = -1,
  RATE_CONSTANT = 0,
  RATE_FILE = 1,
};

/* Special pathway types. */
enum special_pathway_t {
  RFLCT,  /* Special pathway: reflective surface */
  TRANSP, /* Special pathway: transparent surface */
  SINK    /* Special pathway: absorptive surface */
};

struct reaction_def {
  struct sym_entry *sym;
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
  } v;
};

struct reaction_rates {
  struct reaction_rate forward_rate;
  struct reaction_rate backward_rate;
};

MCELL_STATUS
mcell_modify_multiple_rate_constants(struct volume *world, char **names, double *rate_constants, int n_rxns);

MCELL_STATUS
mcell_modify_rate_constant(struct volume *world, char *name, double rate);

MCELL_STATUS
mcell_add_reaction_simplified(
    struct volume *state, 
    struct mcell_species *reactants,
    struct reaction_arrow *arrow,
    struct mcell_species *surfs,
    struct mcell_species *products,
    struct reaction_rates *rates,
    struct sym_entry *pathname);

MCELL_STATUS
mcell_add_reaction(struct notifications *notify,
                   double **r_step_release,
                   struct sym_table_head *rxn_sym_table,
                   u_int radial_subdivisions,
                   double vacancy_search_dist2,
                   struct mcell_species *reactants,
                   struct reaction_arrow *react_arrow,
                   struct mcell_species *surf_class,
                   struct mcell_species *products, struct sym_entry *pathname,
                   struct reaction_rates *rates,
                   const char *forward_rate_filename,
                   const char *backward_rate_filename);

MCELL_STATUS mcell_add_surface_reaction(struct sym_table_head *rxn_sym_table,
                                        int reaction_type,
                                        struct species *surface_class,
                                        struct sym_entry *reactant_sym,
                                        short orient);

MCELL_STATUS
mcell_add_concentration_clamp(struct sym_table_head *rxn_sym_table,
                              struct species *surface_class,
                              struct sym_entry *mol_sym, short orient,
                              double conc);

MCELL_STATUS init_reactions(MCELL_STATE *state);

MCELL_STATUS mcell_change_reaction_rate(MCELL_STATE *state,
                                        const char *reaction_name,
                                        double new_rate);

struct reaction_rates mcell_create_reaction_rates(int forwardRateType,
                                                  double forwardRate,
                                                  int backwardRateType,
                                                  double backwardRate);

struct sym_entry *mcell_new_rxn_pathname(struct volume *state, char *name);
