#ifndef RAN4_H
#define RAN4_H

/*
MCell (tm) Version 2.08 5/18/1998

Copyright (C) 1997,1998, The Salk Institute & Cornell University.
MCell was written jointly by T.M. Bartol Jr. & J.R. Stiles,
with design input for Monte Carlo algorithms from E.E. Salpeter.

  Acknowledgements:
    T.J. Sejnowski for development input and support
    (NSF Grant IBN-9603611), and M.M. Salpeter for fostering
    quantitative experimental applications.  Additional support
    from NIH Grant K08NS01776 (J.R. Stiles).

MCell is a scientific simulation software tool distributed freely to
registered beta-test research laboratories, and must be obtained as a
machine-specific executable program from the authors.  Copying and/or
editing of MCell without the authors' consent is expressly forbidden.
MCell is provided as is and the authors assume no responsibility for
simulation results obtained by users.  Any material published from
MCell simulations must acknowledge the authors and granting
institutions.
*/

void ran4_init(unsigned int *e_state);
double ran4(unsigned int *e_state, double *ran_vec, unsigned int n, double range);
unsigned int iran4(unsigned int *e_state, unsigned int *iran_vec, unsigned int n);
double gaussran4(unsigned int *e_state, double *ran_vec, unsigned int n, double mean, double stddev);

double rng_double(unsigned int index);
unsigned int rng_uint(unsigned int index);

#endif
