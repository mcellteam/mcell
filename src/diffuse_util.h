/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#pragma once

double r_func(double s);
double *init_r_step(int radial_subdivisions);
double *init_r_step_surface(int radial_subdivisions);
double *init_r_step_3d_release(int radial_subdivisions);
double *init_d_step(int radial_directions, unsigned int *actual_directions);
