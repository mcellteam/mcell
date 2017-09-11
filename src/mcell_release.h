/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#pragma once

MCELL_STATUS mcell_create_list_release_site(
  MCELL_STATE *state, struct object *parent, char *site_name,
  struct mcell_species *mol, double *x_pos, double *y_pos, double *z_pos, int n_site,
  struct vector3 *diameter, struct object **new_object);

MCELL_STATUS mcell_create_geometrical_release_site(
    MCELL_STATE *state, struct object *parent, char *site_name, int shape,
    struct vector3 *position, struct vector3 *diameter,
    struct mcell_species *mol, double num, int num_type, double release_prob,
    struct release_pattern *rpatp, struct object **new_object);

MCELL_STATUS mcell_start_release_site(MCELL_STATE *state,
                                      struct sym_entry *sym_ptr,
                                      struct object **obj);

MCELL_STATUS mcell_finish_release_site(struct sym_entry *sym_ptr,
                                       struct object **obj);

/* FIXME: some of the functions below should probably not be part of the API
 * but the parser needs them right now */
int set_release_site_concentration(struct release_site_obj *rel_site_obj_ptr,
                                   double conc);

MCELL_STATUS
mcell_create_region_release(MCELL_STATE *state, struct object *parent,
                            struct object *release_on_in, char *site_name,
                            char *reg_name, struct mcell_species *mol,
                            double num, int num_type, double rel_prob,
                            struct release_pattern *rpatp, struct object **new_object);

MCELL_STATUS
mcell_create_region_release_boolean(MCELL_STATE *state, struct object *parent,
                            char *site_name, struct mcell_species *mol,
                            double num, int num_type, double rel_prob,
                            struct release_pattern *rpatp, struct release_evaluator *rel_eval,
                            struct object **new_object);

struct release_pattern *mcell_create_release_pattern(MCELL_STATE *state, char *name, double delay, 
                                 double release_interval, double train_interval,
                                 double train_duration, int number_of_trains);

int mcell_set_release_site_geometry_region(
    MCELL_STATE *state, struct release_site_obj *rel_site_obj_ptr,
    struct object *objp, struct release_evaluator *re);

int check_release_regions(struct release_evaluator *rel, struct object *parent,
                          struct object *instance);

int is_release_site_valid(struct release_site_obj *rel_site_obj_ptr);

struct release_evaluator *
new_release_region_expr_term(struct sym_entry *my_sym);

void set_release_site_constant_number(struct release_site_obj *rel_site_obj_ptr,
                                      double num);

void set_release_site_gaussian_number(struct release_site_obj *rel_site_obj_ptr,
                                      double mean, double stdev);

struct release_evaluator *
new_release_region_expr_binary(struct release_evaluator *reL,
                               struct release_evaluator *reR, int op);

void set_release_site_location(MCELL_STATE *state,
                               struct release_site_obj *rel_site_obj_ptr,
                               struct vector3 *location);

struct sym_entry *existing_region(MCELL_STATE *state,
								  struct sym_entry *obj_symp,
								  char *region_name);
