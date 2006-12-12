#ifndef MCELL_DIFFUSE_UTIL
#define MCELL_DIFFUSE_UTIL

double dgammln(double xx);
double dgser(double aa, double xx);
double dgcf(double aa, double xx);
double dgammp(double aa, double xx);
double derf(double xx);
double inverf(double xx);

double r_func(double s);
double* init_r_step(int radial_subdivisions);
double* init_r_step_surface(int radial_subdivisions);
double* init_r_step_3d_release(int radial_subdivisions);
double* init_d_step(int radial_directions,unsigned int *actual_directions);

#endif
