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

#include <stdio.h>

/* MCell version as a string */
extern char const mcell_version[];

/* Write the credits to a file handle */
void print_credits(FILE *f);

/* Write the version info to a file handle */
void print_version(FILE *f);

/* Write the version info to a file handle */
void print_full_version(FILE *f);
