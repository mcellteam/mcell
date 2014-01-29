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

 In: mpvp: parser state
     sym_ptr: symbol for the release site
     shape: shape for the release site
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
