#ifndef RAN4_H
#define RAN4_H


void ran4_init(unsigned int *e_state);
double ran4(unsigned int *e_state, double *ran_vec, unsigned int n, double range);
unsigned int iran4(unsigned int *e_state, unsigned int *iran_vec, unsigned int n);
double gaussran4(unsigned int *e_state, double *ran_vec, unsigned int n, double mean, double stddev);

double rng_double(unsigned int index);
unsigned int rng_uint(unsigned int index);

#endif
