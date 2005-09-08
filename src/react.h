#ifndef MCELL_REACT
#define MCELL_REACT

#include "mcell_structs.h"

struct rxn* trigger_unimolecular(int hash,struct abstract_molecule *reac);
struct rxn* trigger_bimolecular(int hashA,int hashB,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB);
struct rxn* trigger_intersect(int hashA,struct abstract_molecule *reacA,
  short orientA,struct wall *w);


int test_unimolecular(struct rxn *rx);
double timeof_unimolecular(struct rxn *rx);
int which_unimolecular(struct rxn *rx);
int test_bimolecular(struct rxn *rx,double scaling);
int test_intersect(struct rxn *rx,double scaling);
void check_probs(struct rxn *rx,double t);


int outcome_products(struct wall *w,struct molecule *reac_m,
  struct grid_molecule *reac_g,struct rxn *rx,int path,struct storage *local,
  short orientA,short orientB,double t,struct vector3 *hitpt,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  struct abstract_molecule *moving);
int outcome_unimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reac,double t);
int outcome_bimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB,double t,struct vector3 *hitpt);
int outcome_intersect(struct rxn *rx, int path, struct wall *surface,
  struct abstract_molecule *reac,short orient,double t,struct vector3 *hitpt);

#endif
