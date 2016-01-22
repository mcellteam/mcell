#ifndef INCLUDED_MACROMOLECULE_H
#define INCLUDED_MACROMOLECULE_H

#include "mcell_structs.h"

struct complex_rate {
  struct complex_rate *next; /* link to next rate table */
  char const *name;          /* name of this rate rule table */

  int num_rules;              /* count of rules in this table */
  int num_neighbors;          /* count of clauses in each rule */
  struct species **neighbors; /* species for rate rule clauses */
  int *invert;                /* invert flags for rate rule clauses */
  signed char *orientations;  /* orients for rate rule clauses */
  double *rates;              /* rates for each rule */
};

#endif
