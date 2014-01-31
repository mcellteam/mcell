/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

#ifndef CREATE_RELEASE_SITE_H
#define CREATE_RELEASE_SITE_H
#include "libmcell.h"

struct release_site_obj *new_release_site(MCELL_STATE *state, char *name);
void set_release_site_location(MCELL_STATE *state,
                               struct release_site_obj *rel_site_obj_ptr,
                               struct vector3 *location);
struct object *start_release_site(MCELL_STATE *state,
                                  struct sym_table *sym_ptr);
struct object *finish_release_site(struct sym_table *sym_ptr);
int is_release_site_valid(struct release_site_obj *rel_site_obj_ptr);
int set_release_site_geometry_region(MCELL_STATE *state,
                                     struct release_site_obj *rel_site_obj_ptr,
                                     struct object *objp,
                                     struct release_evaluator *re);
int check_release_regions(struct release_evaluator *rel,
                          struct object *parent,
                          struct object *instance);
/* Set a constant release quantity from this release site, in units of
 * molecules. */
void set_release_site_constant_number(struct release_site_obj *rel_site_obj_ptr,
                                      double num);
/* Set a gaussian-distributed release quantity from this release site, in units
 * of molecules. */
void set_release_site_gaussian_number(
  struct release_site_obj *rel_site_obj_ptr,
  double mean,
  double stdev);
/* Set a release quantity from this release site based on a fixed concentration
 * in a sphere of a gaussian-distributed diameter with a particular mean and
 * std. deviation. */
void set_release_site_volume_dependent_number(
  struct release_site_obj *rel_site_obj_ptr,
  double mean,
  double stdev,
  double conc);
/* Set a release quantity from this release site based on a fixed concentration
 * within the release-site's area. */
int set_release_site_concentration(struct release_site_obj *rel_site_obj_ptr,
                                   double conc);
/* Set a release quantity from this release site based on a fixed density
 * within the release-site's area. */
int set_release_site_density(struct release_site_obj *rel_site_obj_ptr,
                             double dens);

#endif
