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

#include "geometry_utils.h"
#include "geometry_utils.inl"
#include "rxn_utils.inl"
#include "wall_utils.inl"

namespace MCell {
namespace RegionUtil {


static bool is_any_rxn_reflect_or_absorb_region_border(const BNG::RxnClassesVector& rxns) {
  for (const BNG::RxnClass* rxn_class: rxns) {
    if (rxn_class->is_reflect() || rxn_class->is_absorb_region_border()) {
      return true;
    }
  }

  return false;
}


/***********************************************************************
are_restricted_regions_for_species_on_object:
  In: object
      surface molecule
  Out: true - if there are regions that are restrictive (REFL/ABSORB)
       to the surface molecule on this object
       false - if no such regions found
************************************************************************/
static bool are_restricted_regions_for_species_on_object(
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

  return is_any_rxn_reflect_or_absorb_region_border(matching_rxns);
}


/***********************************************************************
find_restricted_regions_by_object:
  In: object
      surface molecule
  Out: an object's region list that are restrictive (REFL/ABSORB)
       to the surface molecule
       NULL - if no such regions found
  Note: regions called "ALL" or the ones that have ALL_ELEMENTS are not
        included in the return "region list".
************************************************************************/
static void find_restricted_regions_by_object(
    Partition& p, const GeometryObject& obj, const Molecule& sm,
    RegionIndicesSet& res) {

  res.clear();

  const BNG::Species& s = p.get_all_species().get(sm.species_id);
  if (!s.can_interact_with_border()) {
    return;
  }

  struct region *rp;
  struct region_list *rlp, *rlps, *rlp_head = NULL;
  int kk, i, wall_idx = INT_MIN;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];


  for (region_index_t ri: obj.regions) {
    if (ri == obj.encompassing_region_index) {
      continue;
    }

    // find any wall that belongs to this region
    const Region& reg = p.get_region(ri);
    if (reg.walls_and_edges.empty()) {
      continue;
    }

    // we care only about reactive surfaces
    if (!reg.has_surface_class()) {
      continue;
    }

    BNG::RxnClassesVector matching_rxns;
    // NOTE: MCell 3 calls also find_unimol_reactions_with_surf_classes
    // however all reactions should be covered by surface_mol_reactions_with_surf_classes,
    // reactions with surf classes should be always bimolecular - molecule + surf class
    RxnUtil::find_surface_mol_reactions_with_surf_classes(p, sm, obj, matching_rxns);

    if (is_any_rxn_reflect_or_absorb_region_border(matching_rxns)) {
      res.insert(ri);
    }
  }
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
uint determine_molecule_region_topology(
    Partition& p,
    const Molecule* reacA,
    const Molecule* reacB,
    const bool is_unimol,
    RegionIndicesSet& rlp_wall_1,
    RegionIndicesSet& rlp_wall_2,
    RegionIndicesSet& rlp_obj_1,
    RegionIndicesSet& rlp_obj_2) {

  assert(reacA != nullptr);
  const Molecule* sm_1 = (reacA->is_surf()) ? reacA : nullptr;
  const Molecule* sm_2 = (reacB != nullptr && reacB->is_surf()) ? reacB : nullptr;

  uint sm_bitmask = 0;
  rlp_wall_1.clear();
  rlp_wall_2.clear();
  rlp_obj_1.clear();
  rlp_obj_2.clear();

  /* bimolecular surf-surf reaction */
  if (sm_1 != NULL && sm_2 != NULL) {
    const BNG::Species& species1 = p.get_all_species().get(sm_1->species_id);
    const BNG::Species& species2 = p.get_all_species().get(sm_2->species_id);

    const Wall& w_1 = p.get_wall(sm_1->s.wall_index);
    const Wall& w_2 = p.get_wall(sm_2->s.wall_index);

    const GeometryObject& go_1 = p.get_geometry_object(w_1.object_index);
    const GeometryObject& go_2 = p.get_geometry_object(w_2.object_index);

    /* both reactants have restrictive region borders */
    if (
        species1.can_interact_with_border() &&
        are_restricted_regions_for_species_on_object(p, go_1, *sm_1) &&
        species2.can_interact_with_border() &&
        are_restricted_regions_for_species_on_object(p, go_2, *sm_2)
     ) {

      WallUtil::find_restricted_regions_by_wall(p, w_1, *sm_1, rlp_wall_1);
      WallUtil::find_restricted_regions_by_wall(p, w_2, *sm_2, rlp_wall_2);

      /* both reactants are inside their respective restricted regions */
      if (!rlp_wall_1.empty() && !rlp_wall_2.empty()) {
        sm_bitmask |= ALL_INSIDE;
      }
      /* both reactants are outside their respective restricted regions */
      else if (rlp_wall_1.empty() && rlp_wall_2.empty()) {
        find_restricted_regions_by_object(p, go_1, *sm_1, rlp_obj_1);
        find_restricted_regions_by_object(p, go_2, *sm_2, rlp_obj_2);
        sm_bitmask |= ALL_OUTSIDE;
      }
      /* grid1 is inside and grid2 is outside of its respective
       * restrictive region */
      else if (!rlp_wall_1.empty() && rlp_wall_2.empty()) {
        find_restricted_regions_by_object(p, go_2, *sm_2, rlp_obj_2);
        sm_bitmask |= SURF1_IN_SURF2_OUT;
      }
      /* grid2 is inside and grid1 is outside of its respective
       * restrictive region */
      else if (rlp_wall_1.empty() && !rlp_wall_2.empty()) {
        find_restricted_regions_by_object(p, go_1, *sm_1, rlp_obj_1);
        sm_bitmask |= SURF1_OUT_SURF2_IN;
      }
      else {
        assert(false);
      }
    }

    /* only reactant sm_1 has restrictive region border property */
    else if (
         (species1.can_interact_with_border() &&
          are_restricted_regions_for_species_on_object(p, go_1, *sm_1)) &&
        !(species2.can_interact_with_border() &&
          are_restricted_regions_for_species_on_object(p, go_2, *sm_2))
    ){
      WallUtil::find_restricted_regions_by_wall(p, w_1, *sm_1, rlp_wall_1);
      if (!rlp_wall_1.empty()) {
        sm_bitmask |= SURF1_IN;
      }
      else {
        find_restricted_regions_by_object(p, go_1, *sm_1, rlp_obj_1);
        sm_bitmask |= SURF1_OUT;
      }
    }

    /* only reactant "sm_2" has restrictive region border property */
    else if (
        !(species1.can_interact_with_border() &&
          are_restricted_regions_for_species_on_object(p, go_1, *sm_1)) &&
         (species2.can_interact_with_border() &&
          are_restricted_regions_for_species_on_object(p, go_2, *sm_2))
    ){
      WallUtil::find_restricted_regions_by_wall(p, w_2, *sm_2, rlp_wall_2);
      if (!rlp_wall_2.empty()) {
        sm_bitmask |= SURF2_IN;
      }
      else {
        find_restricted_regions_by_object(p, go_2, *sm_2, rlp_obj_2);
        sm_bitmask |= SURF2_OUT;
      }
    }
  }

  /* unimolecular reactions */
  else if ((sm_1 != NULL) && is_unimol) {
    const BNG::Species& species1 = p.get_all_species().get(sm_1->species_id);
    const Wall& w_1 = p.get_wall(sm_1->s.wall_index);
    const GeometryObject& go_1 = p.get_geometry_object(w_1.object_index);

    if (
        (species1.can_interact_with_border() &&
         are_restricted_regions_for_species_on_object(p, go_1, *sm_1))
    ){
      WallUtil::find_restricted_regions_by_wall(p, w_1, *sm_1, rlp_wall_1);
      if (!rlp_wall_1.empty()) {
        sm_bitmask |= ALL_INSIDE;
      }
      else {
        find_restricted_regions_by_object(p, go_1, *sm_1, rlp_obj_1);
        sm_bitmask |= ALL_OUTSIDE;
      }
    }
  }

  return sm_bitmask;
}


/***********************************************************************
 *
 * this function tests if wall target can be reached for product placement
 * based on the previously stored reactant topology based on
 * sm_bitmask. Below, wall 1 is the wall containing reactant 1 and
 * wall 2 is the wall containing reactant 2.
 *
 * in: wall to test for product placement
 *     pointer to array with regions that contain wall 1
 *     pointer to array with regions that contain wall 2
 *     pointer to array with regions that do not contain wall 1
 *     pointer to array with regions that do not contain wall 2
 *
 * out: returns true or false depending if wall target can be
 *      used for product placement.
 *
 ***********************************************************************/
bool product_tile_can_be_reached(
    const Partition& p,
    const wall_index_t wall_index,
    const bool is_unimol,
    const uint sm_bitmask,
    const RegionIndicesSet& rlp_wall_1,
    const RegionIndicesSet& rlp_wall_2,
    const RegionIndicesSet& rlp_obj_1,
    const RegionIndicesSet& rlp_obj_2) {

  bool status = true;

  const Wall& target = p.get_wall(wall_index);

  if (sm_bitmask & ALL_INSIDE) {
    if (is_unimol) {
      if (!WallUtil::wall_belongs_to_all_regions_in_region_list(target, rlp_wall_1)) {
        status = false;
      }
    } else {
      /* bimol reaction */
      if (!WallUtil::wall_belongs_to_all_regions_in_region_list(target, rlp_wall_1) ||
          !WallUtil::wall_belongs_to_all_regions_in_region_list(target, rlp_wall_2)) {
        status = false;
      }
    }
  } else if (sm_bitmask & ALL_OUTSIDE) {
    if (is_unimol) {
      if (WallUtil::wall_belongs_to_any_region_in_region_list(target, rlp_obj_1)) {
        status = false;
      }
    } else {
      if (WallUtil::wall_belongs_to_any_region_in_region_list(target, rlp_obj_1) ||
          WallUtil::wall_belongs_to_any_region_in_region_list(target, rlp_obj_2)) {
        status = false;
      }
    }
  } else if (sm_bitmask & SURF1_IN_SURF2_OUT) {
    if (!WallUtil::wall_belongs_to_all_regions_in_region_list(target, rlp_wall_1) ||
        WallUtil::wall_belongs_to_any_region_in_region_list(target, rlp_obj_2)) {
      status = false;
    }
  } else if (sm_bitmask & SURF1_OUT_SURF2_IN) {
    if (WallUtil::wall_belongs_to_any_region_in_region_list(target, rlp_obj_1) ||
        !WallUtil::wall_belongs_to_all_regions_in_region_list(target, rlp_wall_2)) {
      status = false;
    }
  } else if (sm_bitmask & SURF1_IN) {
    if (!WallUtil::wall_belongs_to_all_regions_in_region_list(target, rlp_wall_1)) {
      status = false;
    }
  } else if (sm_bitmask & SURF1_OUT) {
    if (WallUtil::wall_belongs_to_any_region_in_region_list(target, rlp_obj_1)) {
      status = false;
    }
  } else if (sm_bitmask & SURF2_IN) {
    if (!WallUtil::wall_belongs_to_all_regions_in_region_list(target, rlp_wall_2)) {
      status = false;
    }
  } else if (sm_bitmask & SURF2_OUT) {
    if (WallUtil::wall_belongs_to_any_region_in_region_list(target, rlp_obj_2)) {
      status = false;
    }
  }

  return status;
}


} // namespace RegionUtil
} // namespace MCell
