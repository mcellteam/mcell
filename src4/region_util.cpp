/******************************************************************************
 *
 * Copyright (C) 2021 by
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


#include "region_util.h"

#include "partition.h"

#include "rxn_utils.inl"

namespace MCell {
namespace RegionUtil {


/***********************************************************************
are_restricted_regions_for_species_on_object:
  In: object
      surface molecule
  Out: true - if there are regions that are restrictive (REFL/ABSORB)
       to the surface molecule on this object
       false - if no such regions found
************************************************************************/
bool are_restricted_regions_for_species_on_object(
    Partition& p,
    const GeometryObject& obj,
    const Molecule& sm) {
  assert(sm.is_surf());

  wall_index_t wall_idx = WALL_INDEX_INVALID;
  const BNG::Species& s = p.get_all_species().get(sm.species_id);

  if (!s.can_interact_with_border()) {
    return false;
  }

  // we must check all regions belonging to this object, not just the wall,
  // and get all applicable reactions
  BNG::RxnClassesVector matching_rxns;
  RxnUtil::find_surface_mol_reactions_with_surf_classes(p, sm, obj, matching_rxns);

#if 0
    if (num_matching_rxns > 0) {
      for (int kk = 0; kk < num_matching_rxns; kk++) {
        if ((matching_rxns[kk]->n_pathways == RX_REFLEC) ||
            (matching_rxns[kk]->n_pathways == RX_ABSORB_REGION_BORDER)) {
          return 1;
        }
      }
    }
#endif
  return false;
}


/************************************************************************
 *
 * this function determines where reactants grid1 and grid2 are located
 * (inside/outside) with respect to their restrictive region borders if
 * they have any.
 *
 * in: surface molecule 1 (located on wall 1)
 *     surface molecule 2 (located on wall 2)
 *     pointer to array with restrictive regions which contain wall 1
 *     pointer to array with restrictive regions which contain wall 2
 *     pointer to array with restrictive regions which don't contain wall 1
 *     pointer to array with restrictive regions which don't contain wall 2
 *
 * out: the 4 arrays with pointers to restrictive regions will be filled
 *      and returned
 *
 ***********************************************************************/
int determine_molecule_region_topology(
    const Partition& p,
    const Molecule* reacA,
    const Molecule* reacB,
    const bool is_unimol,
    RegionIndicesSet& rlp_wall_1,
    RegionIndicesSet& rlp_wall_2,
    RegionIndicesSet& rlp_obj_1,
    RegionIndicesSet& rlp_obj_2) {

/*
    struct volume *world, struct surface_molecule *sm_1,
    struct surface_molecule *sm_2, struct region_list **rlp_wall_1_ptr,
    struct region_list **rlp_wall_2_ptr, struct region_list **rlp_obj_1_ptr,
    struct region_list **rlp_obj_2_ptr, bool is_unimol) {*/

  assert(reacA != nullptr);
  const Molecule* sm_1 = (reacA->is_surf()) ? reacA : nullptr;
  const Molecule* sm_2 = (reacB != nullptr && reacB->is_surf()) ? reacB : nullptr;

  const Wall& w_1 = p.get_wall(sm_1->s.wall_index);
  const Wall& w_2 = p.get_wall(sm_2->s.wall_index);

  const GeometryObject& go_1 = p.get_geometry_object(w_1.object_index);
  const GeometryObject& go_2 = p.get_geometry_object(w_2.object_index);

  int sm_bitmask = 0;
  /*struct wall *w_1, *w_2;
  struct region_list *rlp_head_wall_1 = NULL;
  struct region_list *rlp_head_wall_2 = NULL;
  struct region_list *rlp_head_obj_1 = NULL;
  struct region_list *rlp_head_obj_2 = NULL;*/

#if 0
  /* bimolecular surf-surf reaction */
  if (sm_1 != NULL && sm_2 != NULL) {
    const BNG::Species& species1 = p.get_all_species().get(sm_1->species_id);
    const BNG::Species& species2 = p.get_all_species().get(sm_2->species_id);

    /* both reactants have restrictive region borders */
    if (species1.can_interact_with_border() && species2.can_interact_with_border() &&
        are_restricted_regions_for_species_on_object(
            p, go_1, sm_1) &&
        are_restricted_regions_for_species_on_object(
            world, sm_2->grid->surface->parent_object, sm_2)) {
      w_1 = sm_1->grid->surface;
      w_2 = sm_2->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);

      /* both reactants are inside their respective restricted regions */
      if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL)) {
        sm_bitmask |= ALL_INSIDE;
      }
      /* both reactants are outside their respective restricted regions */
      else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL)) {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_bitmask |= ALL_OUTSIDE;
      }
      /* grid1 is inside and grid2 is outside of its respective
       * restrictive region */
      else if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 == NULL)) {
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_bitmask |= SURF1_IN_SURF2_OUT;
      }
      /* grid2 is inside and grid1 is outside of its respective
       * restrictive region */
      else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL)) {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        sm_bitmask |= SURF1_OUT_SURF2_IN;
      }
    }

    /* only reactant sm_1 has restrictive region border property */
    else if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
             are_restricted_regions_for_species_on_object(
                 world, sm_1->grid->surface->parent_object, sm_1) &&
             (!(sm_2->properties->flags & CAN_REGION_BORDER) ||
              !are_restricted_regions_for_species_on_object(
                   world, sm_2->grid->surface->parent_object, sm_2))) {
      w_1 = sm_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      if (rlp_head_wall_1 != NULL) {
        sm_bitmask |= SURF1_IN;
      } else {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        sm_bitmask |= SURF1_OUT;
      }
    }

    /* only reactant "sm_2" has restrictive region border property */
    else if ((sm_2->properties->flags & CAN_REGION_BORDER) &&
             are_restricted_regions_for_species_on_object(
                 world, sm_2->grid->surface->parent_object, sm_2) &&
             (!(sm_1->properties->flags & CAN_REGION_BORDER) ||
              !are_restricted_regions_for_species_on_object(
                   world, sm_1->grid->surface->parent_object, sm_1))) {
      w_2 = sm_2->grid->surface;
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);
      if (rlp_head_wall_2 != NULL) {
        sm_bitmask |= SURF2_IN;
      } else {
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_bitmask |= SURF2_OUT;
      }
    }
  }

  /* unimolecular reactions */
  else if ((sm_1 != NULL) && is_unimol) {
    if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(
            world, sm_1->grid->surface->parent_object, sm_1)) {
      w_1 = sm_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      if (rlp_head_wall_1 != NULL) {
        sm_bitmask |= ALL_INSIDE;
      } else {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        sm_bitmask |= ALL_OUTSIDE;
      }
    }
  }

  *rlp_wall_1_ptr = rlp_head_wall_1;
  *rlp_wall_2_ptr = rlp_head_wall_2;
  *rlp_obj_1_ptr = rlp_head_obj_1;
  *rlp_obj_2_ptr = rlp_head_obj_2;
#endif
  return sm_bitmask;
}

} // namespace RegionUtil
} // namespace MCell
