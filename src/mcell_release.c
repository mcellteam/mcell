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

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "sym_table.h"
#include "logging.h"
#include "mcell_species.h"
#include "mcell_release.h"
#include "mcell_objects.h"

/* static helper functions */
static struct release_site_obj *new_release_site(
    MCELL_STATE *state, char *name);

static struct release_evaluator *
pack_release_expr(
    struct release_evaluator *rel_eval_L,
    struct release_evaluator *rel_eval_R,
    byte op);


/******************************************************************************
 *
 * mcell_create_list_release_site is the main API function for creating
 * release sites from a list of points (LIST based).
 * INPUT:
 * state = world
 * parent = scene
 * site_name
 * mol = really a list of all the mcell_species, constructed using the
 *    inconvenient mcell_add_to_species_list command
 * x,y,z pos
 * number of sites
 * diameter to search for release
 * release object
 *
 ******************************************************************************/
MCELL_STATUS mcell_create_list_release_site(
  MCELL_STATE *state, struct object *parent, char *site_name, struct
  mcell_species *mol, double *x_pos, double *y_pos, double *z_pos, int n_site,
  struct vector3 *diameter, struct object **new_obj) {

  // Create a qualified object name
  // Note: stolen from below. How does this yield a qualified name?
  char qualified_name[20];
    sprintf(qualified_name, "Scene.%s", site_name);

    // Make a new release object - this is later copied into new_obj
  int error_code = 0;
  struct dyngeom_parse_vars *dg_parse = state->dg_parse;
  struct object *release_object = make_new_object(
      dg_parse, state->obj_sym_table, qualified_name, &error_code);

  // Add it to the scene
  // Set the parent of the object to be the root object. Not reciprocal until
  // add_child_objects is called.
  release_object->parent = parent;
  add_child_objects(parent, release_object, release_object);

  /////
  // "Start" the new release site - just like we do in MDL land
  /////
  struct object *dummy = NULL;
  mcell_start_release_site(state, release_object->sym, &dummy);

  // Get a release_site_obj
  struct release_site_obj *releaser = (struct release_site_obj *)release_object->contents;

    // Set it to be a list shape (DOES THIS MATTER?)
  releaser->release_shape = SHAPE_LIST;

  // Set the initial count
  releaser->release_number = 0;

  // Set the diameter - this is the diameter it uses to search for where to place surface mols
  // Note: can be NULL => search diameter = 0, but this means we likely fail to place a lot of surface mols
  releaser->diameter = CHECKED_MALLOC_STRUCT(struct vector3, "release site diameter");
  if (releaser->diameter == NULL) 
  {
    /*free(qualified_name);*/
    return MCELL_FAIL;
  }
  releaser->diameter->x = diameter->x * state->r_length_unit;
  releaser->diameter->y = diameter->y * state->r_length_unit;
  releaser->diameter->z = diameter->z * state->r_length_unit;

  // Go through all the molecules
  int i_site=0;
  struct release_single_molecule *curr_rsm = NULL;
  struct release_single_molecule *rsm = NULL;
  while(i_site < n_site)
  {
    rsm = CHECKED_MALLOC_STRUCT(struct release_single_molecule, "release site molecule position");
    if (rsm == NULL) {
      // Out of memory reading molecule positions
      return MCELL_FAIL;
    }
    
    // Set release properties
    rsm->orient = mol->orient;
    rsm->loc.x = x_pos[i_site] * state->r_length_unit;
    rsm->loc.y = y_pos[i_site] * state->r_length_unit;
    rsm->loc.z = z_pos[i_site] * state->r_length_unit;
    rsm->mol_type = (struct species *)(mol->mol_type->value);
    rsm->next = NULL;

    // mcell_log_raw("Released: %d %s\n", i_site, mol->mol_type->name);

    // If its the first mol, start off the list
    if (i_site == 0)
    {
      releaser->mol_list = rsm; // The head
      curr_rsm = rsm; // For linking the list
    }
    else
    {
      curr_rsm->next = rsm; // Link the previous element to this one
      curr_rsm = rsm; // Advance by one to the tail of the list in anticipation of more linking
    }

    releaser->release_number++; // The count

    // Go to the next site
    i_site++;

    // Make sure to grab the next mol
    mol = mol->next;
  }

  /////
  // "Finish" the release site - yay?
    /////
  mcell_finish_release_site(release_object->sym, &dummy);

  // Copy the new release object into what was passed into the function
  // There probably is something better to do here?
  *new_obj = release_object;

  return MCELL_SUCCESS;
}

/******************************************************************************
 *
 * mcell_create_geometrical_release_site is the main API function for creating
 * geometrical release sites (SHAPE based).
 *
 ******************************************************************************/
MCELL_STATUS mcell_create_geometrical_release_site(
    MCELL_STATE *state, struct object *parent, char *site_name, int shape,
    struct vector3 *position, struct vector3 *diameter,
    struct mcell_species *mol, double num, int num_type, double rel_prob,
    struct release_pattern *rpatp, struct object **new_obj) {

  assert(shape != SHAPE_REGION && shape != SHAPE_LIST);
  assert((((struct species *)mol->mol_type->value)->flags & NOT_FREE) == 0);

  // create qualified object name
  // ecc replaced by sprintf for swig function (macros are bad)
  char *qualified_name = CHECKED_SPRINTF("%s.%s", parent->sym->name, site_name);
  /*char qualified_name[20];*/
  /*sprintf(qualified_name, "Scene.%s", site_name);*/

  int error_code = 0;
  struct dyngeom_parse_vars *dg_parse = state->dg_parse;
  struct object *release_object = make_new_object(
      dg_parse,
      state->obj_sym_table,
      qualified_name,
      &error_code);
  // release_object->parent = state->root_instance;

  // Set the parent of the object to be the root object. Not reciprocal until
  // add_child_objects is called.
  release_object->parent = parent;
  add_child_objects(parent, release_object, release_object);

  struct object *dummy = NULL;
  mcell_start_release_site(state, release_object->sym, &dummy);

  // release site geometry and locations
  struct release_site_obj *releaser =
      (struct release_site_obj *)release_object->contents;
  releaser->release_shape = shape;
  set_release_site_location(state, releaser, position);

  releaser->diameter =
      CHECKED_MALLOC_STRUCT(struct vector3, "release site diameter");
  if (releaser->diameter == NULL) {
    /*free(qualified_name);*/
    return MCELL_FAIL;
  }
  releaser->diameter->x = diameter->x * state->r_length_unit;
  releaser->diameter->y = diameter->y * state->r_length_unit;
  releaser->diameter->z = diameter->z * state->r_length_unit;

  // release probability and release patterns
  if (rel_prob < 0 || rel_prob > 1) {
    /*free(qualified_name);*/
    return MCELL_FAIL;
  }

  if (rpatp != NULL) {
    releaser->pattern = rpatp;
  } else {
    releaser->release_prob = rel_prob;
  }

  /* molecule and molecule number */
  if (num_type == 0) {
    set_release_site_constant_number(releaser, num);
  } else if (num_type == 1) {
    set_release_site_concentration(releaser, num);
  } else {
    return MCELL_FAIL;
  }

  releaser->mol_type = (struct species *)mol->mol_type->value;
  releaser->orientation = mol->orient;

  mcell_finish_release_site(release_object->sym, &dummy);

  *new_obj = release_object;


  //ecc removed for swig function 
  //  free(qualified_name);
 
  return MCELL_SUCCESS;
}

/**************************************************************************
 start_release_site:
  Start parsing the innards of a release site.

 In: state: system state
   sym_ptr: symbol for the release site
 Out: 0 on success, 1 on failure
**************************************************************************/
MCELL_STATUS mcell_start_release_site(MCELL_STATE *state,
                    struct sym_entry *sym_ptr,
                    struct object **obj) {

  struct object *obj_ptr = (struct object *)sym_ptr->value;
  obj_ptr->object_type = REL_SITE_OBJ;
  obj_ptr->contents = new_release_site(state, sym_ptr->name);
  if (obj_ptr->contents == NULL) {
  return MCELL_FAIL;
  }

  *obj = obj_ptr;

  return MCELL_SUCCESS;
}

/**************************************************************************
 finish_release_site:
  Finish parsing the innards of a release site.

 In: sym_ptr: symbol for the release site
 Out: the object, on success, or NULL on failure
**************************************************************************/
MCELL_STATUS mcell_finish_release_site(struct sym_entry *sym_ptr,
                     struct object **obj) {

  struct object *obj_ptr_new = (struct object *)sym_ptr->value;
  no_printf("Release site %s defined:\n", sym_ptr->name);
  if (is_release_site_valid((struct release_site_obj *)obj_ptr_new->contents)) {
  return MCELL_FAIL;
  }
  *obj = obj_ptr_new;

  return MCELL_SUCCESS;
}

/******************************************************************************
 *
 * mcell_create_region_release is the main API function for creating release
 * sites on regions.
 *
 ******************************************************************************/
MCELL_STATUS
mcell_create_region_release(MCELL_STATE *state, struct object *parent,
              struct object *release_on_in, char *site_name,
              char *reg_name, struct mcell_species *mol,
              double num, int num_type, double rel_prob,
              struct release_pattern *rpatp, struct object **new_obj) {

  // create qualified release object name
  char *qualified_name = CHECKED_SPRINTF("%s.%s", parent->sym->name, site_name);

  int error_code = 0;
  struct dyngeom_parse_vars *dg_parse = state->dg_parse;
  struct object *release_object = make_new_object(
      dg_parse,
      state->obj_sym_table,
      qualified_name,
      &error_code);

  // Set the parent of the object to be the root object. Not reciprocal until
  // add_child_objects is called.
  release_object->parent = parent;
  add_child_objects(parent, release_object, release_object);

  struct object *dummy = NULL;
  mcell_start_release_site(state, release_object->sym, &dummy);

  struct release_site_obj *releaser =
    (struct release_site_obj *)release_object->contents;

  struct sym_entry *sym_ptr =
    existing_region(state, release_on_in->sym, reg_name);
  struct release_evaluator *rel_eval = new_release_region_expr_term(sym_ptr);
  mcell_set_release_site_geometry_region(state, releaser, release_on_in,
                     rel_eval);

  // release probability and release patterns
  if (rel_prob < 0 || rel_prob > 1) {
  free(qualified_name);
  return MCELL_FAIL;
  }

  if (rpatp != NULL) {
  releaser->pattern = rpatp;
  } else {
  releaser->release_prob = rel_prob;
  }

  /* molecule and molecule number */
  if (num_type == 0) {
  set_release_site_constant_number(releaser, num);
  } else if (num_type == 1) {
  set_release_site_concentration(releaser, num);
  } else {
  return MCELL_FAIL;
  }

  releaser->mol_type = (struct species *)mol->mol_type->value;
  releaser->orientation = mol->orient;

  mcell_finish_release_site(release_object->sym, &dummy);

  *new_obj = release_object;
  free(qualified_name);
  return MCELL_SUCCESS;
}

/******************************************************************************
 *
 * mcell_create_region_release_boolean is the main API function for creating release
 * sites on regions with boolean logic.
 *
 ******************************************************************************/
MCELL_STATUS
mcell_create_region_release_boolean(MCELL_STATE *state, struct object *parent,
              char *site_name, struct mcell_species *mol,
              double num, int num_type, double rel_prob,
              struct release_pattern *rpatp, struct release_evaluator *rel_eval,
              struct object **new_obj) {

  // create qualified release object name
  char *qualified_name = CHECKED_SPRINTF("%s.%s", parent->sym->name, site_name);

  int error_code = 0;
  struct dyngeom_parse_vars *dg_parse = state->dg_parse;
  struct object *release_object = make_new_object(
      dg_parse,
      state->obj_sym_table,
      qualified_name,
      &error_code);

  // Set the parent of the object to be the root object. Not reciprocal until
  // add_child_objects is called.
  release_object->parent = parent;
  add_child_objects(parent, release_object, release_object);

  struct object *obj_ptr = NULL;
  mcell_start_release_site(state, release_object->sym, &obj_ptr);

  struct release_site_obj *releaser =
    (struct release_site_obj *)release_object->contents;

  mcell_set_release_site_geometry_region(state, releaser, obj_ptr->contents,
                     rel_eval);

  // release probability and release patterns
  if (rel_prob < 0 || rel_prob > 1) {
  free(qualified_name);
  return MCELL_FAIL;
  }

  if (rpatp != NULL) {
  releaser->pattern = rpatp;
  } else {
  releaser->release_prob = rel_prob;
  }

  /* molecule and molecule number */
  if (num_type == 0) {
  set_release_site_constant_number(releaser, num);
  } else if (num_type == 1) {
  set_release_site_concentration(releaser, num);
  } else {
  return MCELL_FAIL;
  }

  releaser->mol_type = (struct species *)mol->mol_type->value;
  releaser->orientation = mol->orient;

  mcell_finish_release_site(release_object->sym, &obj_ptr);

  *new_obj = release_object;
  free(qualified_name);
  return MCELL_SUCCESS;
}

/******************************************************************************
 *
 * mcell_new_release_pattern is the main API function for creating a new
 * release pattern.
 *
 ******************************************************************************/
struct sym_entry *mcell_new_release_pattern(MCELL_STATE *state, char *name) {
  struct sym_entry *st;
  if (retrieve_sym(name, state->rpat_sym_table) != NULL) {
  // Release pattern already defined
  // TO-DO: Ich verlange anstaendige Meldungen, verdammt noch mal!
  // free(name);
  mcell_log_raw("ERROR: Release pattern already defined.");
  return NULL;
  } else if ((st = store_sym(name, RPAT, state->rpat_sym_table,
               NULL)) == NULL) {
  // Out of memory while creating release pattern
  // TO-DO: Ich verlange anstaendige Meldungen, verdammt noch mal!
  // free(name);
  mcell_log_raw("ERROR: Out of memory while creating release pattern.");
  return NULL;
  }

  // free(name);
  return st;
}

/**************************************************************************
* mcell_create_release_pattern:
*    Create a release pattern
**************************************************************************/
struct release_pattern *mcell_create_release_pattern(MCELL_STATE *state, char *name, double delay, 
                 double release_interval, double train_interval,
                 double train_duration, int number_of_trains) {
  struct sym_entry *rpat_sym = mcell_new_release_pattern(state,name);

  if (rpat_sym == NULL) {
  mcell_log_raw("ERROR: Could not create release pattern.");
  return NULL;
  }

  struct release_pattern *rpatp = (struct release_pattern *)rpat_sym->value;

  if (release_interval/state->time_unit <= 0) {
  // Release interval must be set to a positive number
  // TO-DO: Ich verlange anstaendige Meldungen, verdammt noch mal!
  mcell_log_raw("ERROR: Release interval must be set to a positive number.");
  return NULL;
  }
  if (train_interval/state->time_unit <= 0) {
  // Train interval must be set to a positive number
  // TO-DO: Ich verlange anstaendige Meldungen, verdammt noch mal!
  mcell_log_raw("ERROR: Train interval must be set to a positive number.");
  return NULL;
  }
  if (train_duration/state->time_unit > train_interval/state->time_unit) {
  // Train duration must not be longer than the train interval
  // TO-DO: Ich verlange anstaendige Meldungen, verdammt noch mal!
  mcell_log_raw("ERROR: Train duration must not be longer than the train interval.");
  return NULL;
  }
  if (train_duration/state->time_unit <= 0) {
  // Train duration must be set to a positive number
  // TO-DO: Ich verlange anstaendige Meldungen, verdammt noch mal!
  mcell_log_raw("ERROR: Train duration must be set to a positive number.");
  return NULL;
  }

  /* Copy in release pattern */
  if (distinguishable(delay, 0.0, EPS_C)) {
  rpatp->delay = delay/state->time_unit;
  } 
  else {
  rpatp->delay = 0;
  }
  if (distinguishable(release_interval, FOREVER, EPS_C)) {
  rpatp->release_interval = release_interval/state->time_unit;
  } 
  else {
  rpatp->release_interval = FOREVER;
  }
  if (distinguishable(train_interval, FOREVER, EPS_C)) {
  rpatp->train_interval = train_interval/state->time_unit;
  }
  else {
  rpatp->train_interval = FOREVER;
  }
  if (distinguishable(train_duration, FOREVER, EPS_C)) {
  rpatp->train_duration = train_duration/state->time_unit;
  }
  else {
  rpatp->train_duration = FOREVER;
  }
  rpatp->number_of_trains = number_of_trains;

  no_printf("Release pattern %s defined:\n", rpat_sym->name);
  no_printf("\tdelay = %f\n", rpatp->delay);
  no_printf("\trelease_interval = %f\n", rpatp->release_interval);
  no_printf("\ttrain_interval = %f\n", rpatp->train_interval);
  no_printf("\ttrain_duration = %f\n", rpatp->train_duration);
  no_printf("\tnumber_of_trains = %d\n", rpatp->number_of_trains);
  return rpatp;
}


/*************************************************************************
 In: state: system state
   rel_site_obj_ptr: the release site object to validate
   obj_ptr: the object representing this release site
   rel_eval: the release evaluator representing the region of release
 Out: 0 on success, 1 on failure
**************************************************************************/
int mcell_set_release_site_geometry_region(
  MCELL_STATE *state, struct release_site_obj *rel_site_obj_ptr,
  struct object *obj_ptr, struct release_evaluator *rel_eval) {

  rel_site_obj_ptr->release_shape = SHAPE_REGION;
  state->place_waypoints_flag = 1;

  struct release_region_data *rel_reg_data = CHECKED_MALLOC_STRUCT(
    struct release_region_data, "release site on region");
  if (rel_reg_data == NULL) {
  return 1;
  }

  rel_reg_data->n_walls_included = -1; /* Indicates uninitialized state */
  rel_reg_data->cum_area_list = NULL;
  rel_reg_data->wall_index = NULL;
  rel_reg_data->obj_index = NULL;
  rel_reg_data->n_objects = -1;
  rel_reg_data->owners = NULL;
  rel_reg_data->in_release = NULL;
  rel_reg_data->self = obj_ptr;

  rel_reg_data->expression = rel_eval;

  if (check_release_regions(rel_eval, obj_ptr, state->root_instance)) {
  // Trying to release on a region that the release site cannot see! Try
  // grouping the release site and the corresponding geometry with an OBJECT.
  free(rel_reg_data);
  return 2;
  }

  rel_site_obj_ptr->region_data = rel_reg_data;
  return 0;
}

/**************************************************************************
 set_release_site_location:
  Set the location of a release site.

 In: state: system state
   rel_site_obj_ptr: release site
   location: location for release site
 Out: none
**************************************************************************/
void set_release_site_location(MCELL_STATE *state,
                 struct release_site_obj *rel_site_obj_ptr,
                 struct vector3 *location) {
  rel_site_obj_ptr->location = location;
  rel_site_obj_ptr->location->x *= state->r_length_unit;
  rel_site_obj_ptr->location->y *= state->r_length_unit;
  rel_site_obj_ptr->location->z *= state->r_length_unit;
}

/**************************************************************************
 set_release_site_constant_number:
  Set a constant release quantity from this release site, in units of
  molecules.

 In: rel_site_obj_ptr: the release site
   num:  count of molecules to release
 Out: none.  release site object is updated
**************************************************************************/
void set_release_site_constant_number(struct release_site_obj *rel_site_obj_ptr,
                    double num) {
  rel_site_obj_ptr->release_number_method = CONSTNUM;
  rel_site_obj_ptr->release_number = num;
}

/**************************************************************************
 set_release_site_gaussian_number:
  Set a gaussian-distributed release quantity from this release site, in
  units of molecules.

 In: rel_site_obj_ptr: the release site
   mean: mean value of distribution
   stdev: std. dev. of distribution
 Out: none.  release site object is updated
**************************************************************************/
void set_release_site_gaussian_number(struct release_site_obj *rel_site_obj_ptr,
                    double mean, double stdev) {
  rel_site_obj_ptr->release_number_method = GAUSSNUM;
  rel_site_obj_ptr->release_number = mean;
  rel_site_obj_ptr->standard_deviation = stdev;
}

/**************************************************************************
 new_release_region_expr_binary:
  Set the geometry for a particular release site to be a region expression.

 In: parse_state: parser state
   reL:  release evaluation tree (set operations) for left side of expression
   reR:  release evaluation tree for right side of expression
   op:   flags indicating the operation performed by this node
 Out: the release expression, or NULL if an error occurs
**************************************************************************/
struct release_evaluator *
new_release_region_expr_binary(struct release_evaluator *rel_eval_L,
                 struct release_evaluator *rel_eval_R, int op) {
  return pack_release_expr(rel_eval_L, rel_eval_R, op);
}

/*************************************************************************
 check_release_regions:

 In: state:    system state
   rel_eval: an release evaluator (set operations applied to regions)
   parent:   the object that owns this release evaluator
   instance: the root object that begins the instance tree
 Out: 0 if all regions refer to instanced objects or to a common ancestor of
    the object with the evaluator, meaning that the object can be found. 1 if
    any referred-to region cannot be found.
*************************************************************************/
int check_release_regions(struct release_evaluator *rel_eval,
              struct object *parent, struct object *instance) {
  struct object *obj_ptr;

  if (rel_eval->left != NULL) {
  if (rel_eval->op & REXP_LEFT_REGION) {
    obj_ptr =
      common_ancestor(parent, ((struct region *)rel_eval->left)->parent);
    if (obj_ptr == NULL || (obj_ptr->parent == NULL && obj_ptr != instance)) {
    obj_ptr = common_ancestor(instance,
                  ((struct region *)rel_eval->left)->parent);
    }

    if (obj_ptr == NULL) {
    // Region neither instanced nor grouped with release site
    return 2;
    }
  } else if (check_release_regions(rel_eval->left, parent, instance)) {
    return 1;
  }
  }

  if (rel_eval->right != NULL) {
  if (rel_eval->op & REXP_RIGHT_REGION) {
    obj_ptr =
      common_ancestor(parent, ((struct region *)rel_eval->right)->parent);
    if (obj_ptr == NULL || (obj_ptr->parent == NULL && obj_ptr != instance)) {
    obj_ptr = common_ancestor(instance,
                  ((struct region *)rel_eval->right)->parent);
    }

    if (obj_ptr == NULL) {
    // Region not grouped with release site.
    return 3;
    }
  } else if (check_release_regions(rel_eval->right, parent, instance)) {
    return 1;
  }
  }

  return 0;
}

/**************************************************************************
 is_release_site_valid:
  Validate a release site.

 In: rel_site_obj_ptr: the release site object to validate
 Out: 0 if it is valid, 1 if not
**************************************************************************/
int is_release_site_valid(struct release_site_obj *rel_site_obj_ptr) {
  // Unless it's a list release, user must specify MOL type
  if (rel_site_obj_ptr->release_shape != SHAPE_LIST) {
  // Must specify molecule to release using MOLECULE=molecule_name.
  if (rel_site_obj_ptr->mol_type == NULL) {
    return 2;
  }

  // Make sure it's not a surface class
  if ((rel_site_obj_ptr->mol_type->flags & IS_SURFACE) != 0) {
    return 3;
  }
  }

  /* Check that concentration/density status of release site agrees with
   * volume/grid status of molecule */
  if (rel_site_obj_ptr->release_number_method == CCNNUM) {
  // CONCENTRATION may only be used with molecules that can diffuse in 3D.
  if ((rel_site_obj_ptr->mol_type->flags & NOT_FREE) != 0) {
    return 4;
  }
  } else if (rel_site_obj_ptr->release_number_method == DENSITYNUM) {
  // DENSITY may only be used with molecules that can diffuse in 2D.
  if ((rel_site_obj_ptr->mol_type->flags & NOT_FREE) == 0) {
    return 5;
  }
  }

  /* Molecules can only be removed via a region release */
  if (rel_site_obj_ptr->release_shape != SHAPE_REGION &&
    rel_site_obj_ptr->release_number < 0) {
  return 2;
  }

  /* Unless it's a region release we must have a location */
  if (rel_site_obj_ptr->release_shape != SHAPE_REGION) {
  if (rel_site_obj_ptr->location == NULL) {
    // Release site is missing location.
    if (rel_site_obj_ptr->release_shape != SHAPE_LIST ||
      rel_site_obj_ptr->mol_list == NULL) {
    return 6;
    } else {
    // Give it a default location of (0, 0, 0)
    rel_site_obj_ptr->location =
      CHECKED_MALLOC_STRUCT(struct vector3, "release site location");
    if (rel_site_obj_ptr->location == NULL)
      return 1;
    rel_site_obj_ptr->location->x = 0;
    rel_site_obj_ptr->location->y = 0;
    rel_site_obj_ptr->location->z = 0;
    }
  }
  no_printf("\tLocation = [%f,%f,%f]\n", rel_site_obj_ptr->location->x,
        rel_site_obj_ptr->location->y, rel_site_obj_ptr->location->z);
  }
  return 0;
}

/**************************************************************************
 set_release_site_concentration:
  Set a release quantity from this release site based on a fixed
  concentration within the release-site's area.

 In: rel_site_obj_ptr: the release site
   conc: concentration for release
 Out: 0 on success, 1 on failure.  release site object is updated
**************************************************************************/
int set_release_site_concentration(struct release_site_obj *rel_site_obj_ptr,
                   double conc) {
  if (rel_site_obj_ptr->release_shape == SHAPE_SPHERICAL_SHELL) {
  return 1;
  }
  rel_site_obj_ptr->release_number_method = CCNNUM;
  rel_site_obj_ptr->concentration = conc;
  return 0;
}

/**************************************************************************
 mdl_new_release_region_expr_term:
  Create a new "release on region" expression term.

 In: my_sym: the symbol for the region comprising this term in the expression
 Out: the release evaluator on success, or NULL if allocation fails
**************************************************************************/
struct release_evaluator *
new_release_region_expr_term(struct sym_entry *my_sym) {

  struct release_evaluator *rel_eval =
    CHECKED_MALLOC_STRUCT(struct release_evaluator, "release site on region");
  if (rel_eval == NULL) {
  return NULL;
  }

  rel_eval->op = REXP_NO_OP | REXP_LEFT_REGION;
  rel_eval->left = my_sym->value;
  rel_eval->right = NULL;

  ((struct region *)rel_eval->left)->flags |= COUNT_CONTENTS;
  return rel_eval;
}

/******************************************************************************
 *
 * static helper functions
 *
 *****************************************************************************/

/*************************************************************************
 existing_region:
  Find an existing region.  Print an error message if it isn't found.

 In:  obj_symp: object on which to find the region
    name: region name
 Out: the region, or NULL if not found
 NOTE: This is similar to mdl_existing_region
*************************************************************************/
struct sym_entry *existing_region(MCELL_STATE *state,
                  struct sym_entry *obj_symp,
                  char *region_name) {
  char *full_name = CHECKED_SPRINTF("%s,%s", obj_symp->name, region_name);
  if (full_name == NULL) {
  // free(full_name);
  return NULL;
  }

  struct sym_entry *symp = retrieve_sym(full_name, state->reg_sym_table);

  free(full_name);
  return symp;
}

/**************************************************************************
 new_release_site:
  Create a new release site.

 In: state: system state
   name: name for the new site
 Out: an empty release site, or NULL if allocation failed
**************************************************************************/
struct release_site_obj *new_release_site(MCELL_STATE *state, char *name) {
  struct release_site_obj *rel_site_obj_ptr;
  if ((rel_site_obj_ptr = CHECKED_MALLOC_STRUCT(struct release_site_obj,
                        "release site")) == NULL)
  return NULL;
  rel_site_obj_ptr->location = NULL;
  rel_site_obj_ptr->mol_type = NULL;
  rel_site_obj_ptr->release_number_method = CONSTNUM;
  rel_site_obj_ptr->release_shape = SHAPE_UNDEFINED;
  rel_site_obj_ptr->orientation = 0;
  rel_site_obj_ptr->release_number = 0;
  rel_site_obj_ptr->mean_diameter = 0;
  rel_site_obj_ptr->concentration = 0;
  rel_site_obj_ptr->standard_deviation = 0;
  rel_site_obj_ptr->diameter = NULL;
  rel_site_obj_ptr->region_data = NULL;
  rel_site_obj_ptr->mol_list = NULL;
  rel_site_obj_ptr->release_prob = 1.0;
  rel_site_obj_ptr->pattern = state->default_release_pattern;
  struct periodic_image *periodic_box = CHECKED_MALLOC_STRUCT(
  struct periodic_image, "periodic image descriptor");
  rel_site_obj_ptr->periodic_box = periodic_box;
  rel_site_obj_ptr->periodic_box->x = 0;
  rel_site_obj_ptr->periodic_box->y = 0;
  rel_site_obj_ptr->periodic_box->z = 0;
  // if ((rel_site_obj_ptr->name = mdl_strdup(name)) == NULL)
  if ((rel_site_obj_ptr->name = strdup(name)) == NULL) {
  free(rel_site_obj_ptr);
  return NULL;
  }
  return rel_site_obj_ptr;
}

/*************************************************************************
 pack_release_expr:

 In: rel_eval_L:  release evaluation tree (set operations) for left side of
expression
   rel_eval_R:  release evaluation tree for right side of expression
   op:   flags indicating the operation performed by this node
 Out: release evaluation tree containing the two subtrees and the
    operation
 Note: singleton elements (with REXP_NO_OP operation) are compacted by
     this function and held simply as the corresponding region, not
     the NO_OP operation of that region (the operation is needed for
     efficient parsing)
*************************************************************************/
struct release_evaluator *
pack_release_expr(struct release_evaluator *rel_eval_L,
          struct release_evaluator *rel_eval_R, byte op) {

  struct release_evaluator *rel_eval = NULL;

  if ((rel_eval_R->op & REXP_MASK) == REXP_NO_OP &&
    (rel_eval_R->op & REXP_LEFT_REGION) != 0) {
  if ((rel_eval_L->op & REXP_MASK) == REXP_NO_OP &&
    (rel_eval_L->op & REXP_LEFT_REGION) != 0) {
    rel_eval = rel_eval_L;
    rel_eval->right = rel_eval_R->left;
    rel_eval->op = op | REXP_LEFT_REGION | REXP_RIGHT_REGION;
    free(rel_eval_R);
  } else {
    rel_eval = rel_eval_R;
    rel_eval->right = rel_eval->left;
    rel_eval->left = (void *)rel_eval_L;
    rel_eval->op = op | REXP_RIGHT_REGION;
  }
  } else if ((rel_eval_L->op & REXP_MASK) == REXP_NO_OP &&
       (rel_eval_L->op & REXP_LEFT_REGION) != 0) {
  rel_eval = rel_eval_L;
  rel_eval->right = (void *)rel_eval_R;
  rel_eval->op = op | REXP_LEFT_REGION;
  } else {
  rel_eval = CHECKED_MALLOC_STRUCT(struct release_evaluator,
                   "release region expression");
  if (rel_eval == NULL) {
    return NULL;
  }

  rel_eval->left = (void *)rel_eval_L;
  rel_eval->right = (void *)rel_eval_R;
  rel_eval->op = op;
  }

  return rel_eval;
}
