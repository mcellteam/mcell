/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by *
 * The Salk Institute for Biological Studies and *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or *
 * modify it under the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 *
 * of the License, or (at your option) any later version. *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *
 * GNU General Public License for more details. *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License *
 * along with this program; if not, write to the Free Software *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 *USA. *
 *                                                                                 *
 ***********************************************************************************/

#include "create_release_site.h"
#include "create_object.h"
#include "logging.h"
#include "sym_table.h"

#include <stdlib.h>
#include <string.h>

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

  if (!(op & REXP_INCLUSION) && (rel_eval_R->op & REXP_MASK) == REXP_NO_OP &&
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
  } else if (!(op & REXP_INCLUSION) &&
             (rel_eval_L->op & REXP_MASK) == REXP_NO_OP &&
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