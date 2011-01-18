#ifndef MCELL_REACT
#define MCELL_REACT

#include "mcell_structs.h"


/* In react_trig.c */
struct rxn* trigger_unimolecular(u_int hash,struct abstract_molecule *reac);
int trigger_surface_unimol(struct abstract_molecule *reac,struct wall *w, struct rxn **matching_rxns);
int trigger_bimolecular_preliminary(u_int hashA,u_int hashB,
  struct species *reacA,struct species *reacB);
int trigger_trimolecular_preliminary(u_int hashA, u_int hashB, u_int hashC,
  struct species *reacA, struct species *reacB, struct species *reacC);
int trigger_bimolecular(u_int hashA,u_int hashB,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB, struct rxn **matching_rxns);
int trigger_trimolecular(u_int hashA,u_int hashB, u_int hashC,
  struct species *reacA,struct species *reacB,
  struct species *reacC, int orientA, int orientB, int orientC, 
  struct rxn **matching_rxns);
int trigger_intersect(u_int hashA,struct abstract_molecule *reacA,
  short orientA,struct wall *w, struct rxn **matching_rxns, int allow_rx_transp,  int allow_rx_reflec, int allow_rx_absorb_reg_border);


/* In react_cond.c */
double timeof_unimolecular(struct rxn *rx, struct abstract_molecule *a);
double timeof_special_unimol(struct rxn *rxuni,struct rxn *rxsurf, struct abstract_molecule *a);
int which_unimolecular(struct rxn *rx, struct abstract_molecule *a);
int is_surface_unimol(struct rxn *rxuni,struct rxn *rxsurf, struct abstract_molecule *a);
int test_bimolecular(struct rxn *rx, double scaling, double local_prob_factor, struct abstract_molecule *a1, struct abstract_molecule *a2);
int test_many_bimolecular(struct rxn **rx, double *scaling, int n, int *chosen_pathway, struct abstract_molecule **complexes, int *complex_limits);
int test_many_bimolecular_all_neighbors(struct rxn **rx, double *scaling, double local_prob_factor, int n, int *chosen_pathway, struct abstract_molecule **complexes, int *complex_limits);
int test_many_reactions_all_neighbors(struct rxn **rx, double *scaling, double *local_prob_factor, int n, int *chosen_pathway);
int test_intersect(struct rxn *rx,double scaling);
int test_many_intersect(struct rxn **rx,double scaling, int n, int *chosen_pathway);
struct rxn * test_many_intersect_unimol(struct rxn **rx, int n);
void check_probs(struct rxn *rx,double t);


/* In react_outc.c */
int outcome_unimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reac,double t);
int outcome_bimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB,double t,struct vector3 *hitpt,
  struct vector3 *loc_okay);
int outcome_trimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  struct abstract_molecule *reacC, short orientA, short orientB, 
  short orientC, double t,struct vector3 *hitpt,struct vector3 *loc_okay);
int outcome_intersect(struct rxn *rx, int path, struct wall *surface,
  struct abstract_molecule *reac,short orient,double t,struct vector3 *hitpt,
  struct vector3 *loc_okay);
int is_compatible_surface(void *req_species, struct wall *w);

#endif
