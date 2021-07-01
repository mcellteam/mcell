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

// cannot be compiled on MacOS
//MCELL_STATUS
//mcell_modify_multiple_rate_constants(struct volume *world, char **names, double *rate_constants, int n_rxns);


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
mcell_add_clamp(struct sym_table_head *rxn_sym_table,
                              struct species *surface_class,
                              struct sym_entry *mol_sym, short orient,
                              int clamp_type,
                              double clamp_value);                              

MCELL_STATUS init_reactions(MCELL_STATE *state);

MCELL_STATUS mcell_change_reaction_rate(MCELL_STATE *state,
                                        const char *reaction_name,
                                        double new_rate);

struct reaction_rates mcell_create_reaction_rates(int forwardRateType,
                                                  double forwardRate,
                                                  int backwardRateType,
                                                  double backwardRate);

struct sym_entry *mcell_new_rxn_pathname(struct volume *state, char *name);
