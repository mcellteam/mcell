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

MCELL_STATUS mcell_create_list_release_site(
  MCELL_STATE *state, struct geom_object *parent, char *site_name,
  struct mcell_species *mol, double *x_pos, double *y_pos, double *z_pos, int n_site,
  struct vector3 *diameter, struct geom_object **new_object);

MCELL_STATUS mcell_create_geometrical_release_site(
    MCELL_STATE *state, struct geom_object *parent, const char *site_name, int shape,
    struct vector3 *position, struct vector3 *diameter,
    struct mcell_species *mol, double num, int num_type, double release_prob,
    struct release_pattern *rpatp, struct geom_object **new_object);

MCELL_STATUS mcell_start_release_site(MCELL_STATE *state,
                                      struct sym_entry *sym_ptr,
                                      struct geom_object **obj);

MCELL_STATUS mcell_finish_release_site(struct sym_entry *sym_ptr,
                                       struct geom_object **obj);

/* FIXME: some of the functions below should probably not be part of the API
 * but the parser needs them right now */
int set_release_site_concentration(struct release_site_obj *rel_site_obj_ptr,
                                   double conc);

MCELL_STATUS
mcell_create_region_release(MCELL_STATE *state, struct geom_object *parent,
                            struct geom_object *release_on_in, char *site_name,
                            char *reg_name, struct mcell_species *mol,
                            double num, int num_type, double rel_prob,
                            struct release_pattern *rpatp, struct geom_object **new_object);

MCELL_STATUS
mcell_create_region_release_boolean(MCELL_STATE *state, struct geom_object *parent,
                            char *site_name, struct mcell_species *mol,
                            double num, int num_type, double rel_prob,
                            struct release_pattern *rpatp, struct release_evaluator *rel_eval,
                            struct geom_object **new_object);

struct release_pattern *mcell_create_release_pattern(MCELL_STATE *state, char *name, double delay, 
                                 double release_interval, double train_interval,
                                 double train_duration, int number_of_trains);

int mcell_set_release_site_geometry_region(
    MCELL_STATE *state, struct release_site_obj *rel_site_obj_ptr,
    struct geom_object *objp, struct release_evaluator *re);

int check_release_regions(struct release_evaluator *rel, struct geom_object *parent,
                          struct geom_object *instance);

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
