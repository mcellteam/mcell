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

#include "create_release_site.h"
#include "create_object.h"
#include "logging.h"

#include <stdlib.h>
#include <string.h>



/**************************************************************************
 new_release_site:
    Create a new release site.

 In: state: system state
     name: name for the new site
 Out: an empty release site, or NULL if allocation failed
**************************************************************************/
struct release_site_obj *
new_release_site(MCELL_STATE *state, char *name)
{
  struct release_site_obj *rel_site_obj_ptr;
  if ((rel_site_obj_ptr = CHECKED_MALLOC_STRUCT(struct release_site_obj, "release site")) == NULL)
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
  //if ((rel_site_obj_ptr->name = mdl_strdup(name)) == NULL)
  if ((rel_site_obj_ptr->name = strdup(name)) == NULL)
  {
    free(rel_site_obj_ptr);
    return NULL;
  }
  return rel_site_obj_ptr;
}



/**************************************************************************
 set_release_site_location:
    Set the location of a release site.

 In: state: system state
     rel_site_obj_ptr: release site
     location: location for release site
 Out: none
**************************************************************************/
void
set_release_site_location(MCELL_STATE *state,
                          struct release_site_obj *rel_site_obj_ptr,
                          struct vector3 *location)
{
  rel_site_obj_ptr->location = location;
  rel_site_obj_ptr->location->x *= state->r_length_unit;
  rel_site_obj_ptr->location->y *= state->r_length_unit;
  rel_site_obj_ptr->location->z *= state->r_length_unit;
}



/**************************************************************************
 start_release_site:
    Start parsing the innards of a release site.

 In: state: system state
     sym_ptr: symbol for the release site
 Out: 0 on success, 1 on failure
**************************************************************************/
struct object *
start_release_site(MCELL_STATE *state,
                   struct sym_table *sym_ptr)
{
  struct object *obj_ptr = (struct object *) sym_ptr->value;
  obj_ptr->object_type = REL_SITE_OBJ;
  obj_ptr->contents = new_release_site(state, sym_ptr->name);
  if (obj_ptr->contents == NULL) {
    return NULL;
  }

  return obj_ptr;
}



/**************************************************************************
 finish_release_site:
    Finish parsing the innards of a release site.

 In: sym_ptr: symbol for the release site
 Out: the object, on success, or NULL on failure
**************************************************************************/
struct object *
finish_release_site(struct sym_table *sym_ptr)
{
  struct object *obj_ptr_new = (struct object *) sym_ptr->value;
  no_printf("Release site %s defined:\n", sym_ptr->name);
  if (is_release_site_valid((struct release_site_obj *) obj_ptr_new->contents)) {
    return NULL;
  }
  return obj_ptr_new;
}



/**************************************************************************
 is_release_site_valid:
    Validate a release site.

 In: rel_site_obj_ptr: the release site object to validate
 Out: 0 if it is valid, 1 if not
**************************************************************************/
int
is_release_site_valid(struct release_site_obj *rel_site_obj_ptr)
{
  // Unless it's a list release, user must specify MOL type 
  if (rel_site_obj_ptr->release_shape != SHAPE_LIST)
  {
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
  if (rel_site_obj_ptr->release_number_method == CCNNUM)
  {
    // CONCENTRATION may only be used with molecules that can diffuse in 3D.
    if ((rel_site_obj_ptr->mol_type->flags & NOT_FREE) != 0) {
      return 4;
    }
  }
  else if (rel_site_obj_ptr->release_number_method == DENSITYNUM)
  {
    // DENSITY may only be used with molecules that can diffuse in 2D.
    if ((rel_site_obj_ptr->mol_type->flags & NOT_FREE) == 0) {
      return 5;
    }
  }

  /* Unless it's a region release we must have a location */
  if (rel_site_obj_ptr->release_shape != SHAPE_REGION)
  {
    if (rel_site_obj_ptr->location == NULL)
    {
      // Release site is missing location.
      if (rel_site_obj_ptr->release_shape!=SHAPE_LIST || rel_site_obj_ptr->mol_list==NULL) {
        return 6;
      }
      else
      {
        // Give it a default location of (0, 0, 0)
        rel_site_obj_ptr->location = CHECKED_MALLOC_STRUCT(struct vector3, "release site location");
        if (rel_site_obj_ptr->location==NULL)
          return 1;
        rel_site_obj_ptr->location->x = 0;
        rel_site_obj_ptr->location->y = 0;
        rel_site_obj_ptr->location->z = 0;
      }
    }
    no_printf(
      "\tLocation = [%f,%f,%f]\n",
      rel_site_obj_ptr->location->x,
      rel_site_obj_ptr->location->y,
      rel_site_obj_ptr->location->z);
  }
  return 0;
}



/**************************************************************************
 set_release_site_geometry_region:
    Set the geometry for a particular release site to be a region expression.

 In: state: system state
     rel_site_obj_ptr: the release site object to validate
     obj_ptr: the object representing this release site
     rel_eval: the release evaluator representing the region of release
 Out: 0 on success, 1 on failure
**************************************************************************/
int
set_release_site_geometry_region(MCELL_STATE *state,
                                 struct release_site_obj *rel_site_obj_ptr,
                                 struct object *obj_ptr,
                                 struct release_evaluator *rel_eval)
{

  rel_site_obj_ptr->release_shape = SHAPE_REGION;
  state->place_waypoints_flag = 1;

  struct release_region_data *rel_reg_data = CHECKED_MALLOC_STRUCT(
    struct release_region_data,
    "release site on region");
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

  if (check_release_regions(rel_eval, obj_ptr, state->root_instance))
  {
    // Trying to release on a region that the release site cannot see! Try
    // grouping the release site and the corresponding geometry with an OBJECT.
    free(rel_reg_data);
    return 2;
  }

  rel_site_obj_ptr->region_data = rel_reg_data;
  return 0;
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
int
check_release_regions(struct release_evaluator *rel_eval,
                      struct object *parent,
                      struct object *instance)
{
  struct object *obj_ptr;

  if (rel_eval->left != NULL)
  {
    if (rel_eval->op & REXP_LEFT_REGION)
    {
      obj_ptr = common_ancestor(parent, ((struct region*) rel_eval->left)->parent);
      if (obj_ptr == NULL || (obj_ptr->parent == NULL && obj_ptr!=instance)) {
        obj_ptr = common_ancestor(instance, ((struct region*) rel_eval->left)->parent);
      }

      if (obj_ptr == NULL)
      {
        // Region neither instanced nor grouped with release site
        return 2;
      }
    }
    else if (check_release_regions(rel_eval->left, parent, instance)) {
      return 1;
    }
  }

  if (rel_eval->right != NULL)
  {
    if (rel_eval->op & REXP_RIGHT_REGION)
    {
      obj_ptr = common_ancestor(parent, ((struct region*)rel_eval->right)->parent);
      if (obj_ptr == NULL || (obj_ptr->parent == NULL && obj_ptr != instance)) {
        obj_ptr = common_ancestor(instance, ((struct region*)rel_eval->right)->parent);
      }

      if (obj_ptr == NULL)
      {
        // Region not grouped with release site.
        return 3;
      }
    }
    else if (check_release_regions(rel_eval->right, parent, instance)) {
      return 1;
    }
  }

  return 0;
}



/**************************************************************************
 set_release_site_constant_number:
    Set a constant release quantity from this release site, in units of
    molecules.

 In: rel_site_obj_ptr: the release site
     num:  count of molecules to release
 Out: none.  release site object is updated
**************************************************************************/
void 
set_release_site_constant_number(struct release_site_obj *rel_site_obj_ptr,
                                 double num)
{
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
void 
set_release_site_gaussian_number(struct release_site_obj *rel_site_obj_ptr,
                                 double mean,
                                 double stdev)
{
  rel_site_obj_ptr->release_number_method = GAUSSNUM;
  rel_site_obj_ptr->release_number = mean;
  rel_site_obj_ptr->standard_deviation = stdev;
}



/**************************************************************************
 set_release_site_volume_dependent_number:
    Set a release quantity from this release site based on a fixed
    concentration in a sphere of a gaussian-distributed diameter with a
    particular mean and std. deviation.

 In: rel_site_obj_ptr: the release site
     mean: mean value of distribution of diameters
     stdev: std. dev. of distribution of diameters
     conc: concentration for release
 Out: none.  release site object is updated
**************************************************************************/
void 
set_release_site_volume_dependent_number(struct release_site_obj *rel_site_obj_ptr,
                                         double mean,
                                         double stdev,
                                         double conc)
{
  rel_site_obj_ptr->release_number_method = VOLNUM;
  rel_site_obj_ptr->mean_diameter = mean;
  rel_site_obj_ptr->standard_deviation = stdev;
  rel_site_obj_ptr->concentration = conc;
}



/**************************************************************************
 set_release_site_concentration:
    Set a release quantity from this release site based on a fixed
    concentration within the release-site's area.

 In: rel_site_obj_ptr: the release site
     conc: concentration for release
 Out: 0 on success, 1 on failure.  release site object is updated
**************************************************************************/
int 
set_release_site_concentration(struct release_site_obj *rel_site_obj_ptr,
                               double conc)
{
  if (rel_site_obj_ptr->release_shape == SHAPE_SPHERICAL_SHELL) {
    return 1;
  }
  rel_site_obj_ptr->release_number_method = CCNNUM;
  rel_site_obj_ptr->concentration = conc;
  return 0;
}



/**************************************************************************
 set_release_site_density:
    Set a release quantity from this release site based on a fixed
    density within the release-site's area.  (Hopefully we're talking about a
    surface release here.)

 In: rel_site_obj_ptr: the release site
     dens: density for release
 Out: 0 on success, 1 on failure.  release site object is updated
**************************************************************************/
int 
set_release_site_density(struct release_site_obj *rel_site_obj_ptr,
                         double dens)
{

  rel_site_obj_ptr->release_number_method = DENSITYNUM;
  rel_site_obj_ptr->concentration = dens;
  return 0;
}
