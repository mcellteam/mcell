/*
 * rxn_compartment_utils.cpp
 *
 *  Created on: Oct 27, 2020
 *      Author: ahusar
 */

// Code below implements semantic checks defined here:
// https://docs.google.com/document/d/1E5T0dLVXgYLml7VlfoXwv_gAV4T6q29fbOoKu6gXnHo/edit?usp=sharing

#include <vector>

#include "bng/rxn_rule.h"
#include "bng/bng_data.h"
#include "bng/cplx.h"

#include "bng/rxn_compartment_utils.h"

using namespace std;

namespace BNG {

#define CHECK(call) { string res = call; if (res != "") return res; }


static string check_no_in_out_compartment(const RxnRule& r, const CplxVector& substances) {
  // 5. It is not allowed to use @IN or @OUT in a volume reaction
  for (const Cplx& subst: substances) {
    if (subst.has_compartment_class_in_out()) {
      return "Reaction rule " + r.to_str(false, false, false) + ": compartment class " +
          compartment_id_to_str(subst.get_primary_compartment_id()) + " is not allowed in volume reactions.";
    }
  }
  return "";
}



static bool has_vol_substance(const CplxVector& substances) {
  // 5. It is not allowed to use @IN or @OUT in a volume reaction
  for (const Cplx& subst: substances) {
    if (subst.is_vol()) {
      return true;
    }
  }
  return false;
}


static bool has_only_simple_substances(const CplxVector& substances) {
  for (const Cplx& subst: substances) {
    if (!subst.is_simple()) {
      return false;
    }
  }
  return true;
}


static bool all_have_compartment(const CplxVector& substances) {
  // 5. It is not allowed to use @IN or @OUT in a volume reaction
  for (const Cplx& subst: substances) {
    if (!subst.has_compartment()) {
      return false;
    }
  }
  return true;
}


static string check_surf_have_no_compartment(const RxnRule& r, const CplxVector& substances) {
  // 5. It is not allowed to use @IN or @OUT in a volume reaction
  for (const Cplx& subst: substances) {
    if (subst.is_surf() && subst.has_compartment()) {
      return "In a surface reaction rule " + r.to_str(false, false, false) + " that uses compartments " +
          COMPARTMENT_NAME_IN + " or " + COMPARTMENT_NAME_OUT + ", no surface products or reactants must have their " +
          "compartment specified, error for " + subst.to_str() + ".";
    }
  }
  return "";
}


static string check_vol_have_in_out_compartment(const RxnRule& r, const CplxVector& substances) {
  // 5. It is not allowed to use @IN or @OUT in a volume reaction
  for (const Cplx& subst: substances) {
    if (subst.is_vol() && !subst.has_compartment_class_in_out()) {
      return "In a surface reaction rule " + r.to_str(false, false, false) + " that uses compartments " +
          COMPARTMENT_NAME_IN + " or " + COMPARTMENT_NAME_OUT + ", no surface products or reactants must have their " +
          "compartment specified, error for " + subst.to_str() + ".";
    }
  }
  return "";
}


// surf_comp is set to a value other than COMPARTMENT_ID_NONE when the rxn
// has surface reactants that use specific compartment,
// also sets uses_in_out_compartment to true if the rxn uses @IN or @OUT
// (only surface rxns may use these compartment classes)
static string check_surface_compartments(
    const BNGData& bng_data, const RxnRule& r,
    compartment_id_t& surf_comp_id, bool& uses_in_out_compartment, bool& has_surf_mols_and_vol_compartments) {
  // NOTE: this seems to belongs to the BNGL library but we do not know
  // whether the molecules types are surface or volume there

  surf_comp_id = COMPARTMENT_ID_NONE;
  uses_in_out_compartment = false;
  has_surf_mols_and_vol_compartments = false;

  // this is applicable only when there is a surface reactant
  if (!r.is_surf_rxn()) {
    // 4. A volume reaction must not have surface products
    for (const Cplx& prod: r.products) {
      if (prod.is_surf()) {
        return "Reaction rule " + r.to_str(false, false, false) + " has no surface reactants but has surface "
            "products, this is not allowed because it is not known where to create these surface products.";
      }
    }
    // no IN/OUT compartments are used
    CHECK(check_no_in_out_compartment(r, r.reactants));
    CHECK(check_no_in_out_compartment(r, r.products));

    // nothing to do otherwise
    return "";
  }

  // do we have a surface compartment specified?
  bool has_surf_reactants = false;
  bool uses_vol_compartments = false;
  vector<compartment_id_t> surf_compartments;
  for (size_t i = 0; i < r.reactants.size(); i++) {
    const Cplx& reac = r.reactants[i];

    if (reac.is_surf()) {
      has_surf_reactants = true;
      // remember the first surface compartment
      surf_compartments.push_back(reac.get_primary_compartment_id());
    }

    uint_set<compartment_id_t> used_compartments;
    reac.get_used_compartments(used_compartments);
    for (compartment_id_t cid: used_compartments) {

      if (is_specific_compartment_id(cid) && bng_data.get_compartment(cid).is_3d) {
        uses_vol_compartments = true;
      }
    }

    uses_in_out_compartment = uses_in_out_compartment || reac.has_compartment_class_in_out();
  }

  has_surf_mols_and_vol_compartments = has_surf_reactants && uses_vol_compartments;

  if (!uses_in_out_compartment) {
    // 3. In a surface reaction:
    //  a. if reaction does not use @IN/@OUT and a volume reactant or product is present:
    //     all substances must have their specific compartments specified
    //  b. if reaction does not use @IN/@OUT:
    //     all surface substances have the same compartment or no compartment at all

    //  c. check for cases that we cannot handle yet because we have a single compartment for
    //     a complex:
    //
    // S + * -> V + *
    //    - allowed only:
    //        - when all substances are simple - we cannot distinguish
    //          it from MCell reaction S; + * -> V; + *
    //          or
    //        - when compartments are specified for all substances
    //
    bool has_vol_product = has_vol_substance(r.products);

    if (has_vol_product) {

      bool only_simple = has_only_simple_substances(r.reactants) && has_only_simple_substances(r.products);
      bool all_compartments = all_have_compartment(r.reactants) && all_have_compartment(r.reactants);

      if (!(only_simple || all_compartments)) {
        return "Cannot determine orientation or compartment of a volume product of a surface reaction rule " + r.to_str(false, false, false) + ". "
            "Legal BNGL reactions in the form S(s!1).V(v!1) -> V(s) + S(s) are not supported yet. " +
            "Either specify compartments for each reactant and product or use compartment class "
            "@" + COMPARTMENT_NAME_IN + " or @" + COMPARTMENT_NAME_OUT + " for the volume products "
            "(these compartment classes are however not supported by official BioNetGen tool yet).";
      }
    }

    assert(surf_compartments.size() <= 2);
    if (surf_compartments.size() == 1) {
      surf_comp_id = surf_compartments[0];
    }
    else if (surf_compartments.size() == 2 && !r.is_intermembrane_surf_rxn()) {
      if (surf_compartments[0] != surf_compartments[1]) {
        return "Reaction rule " + r.to_str(false, false, false) + ": all compartments of surface molecules on the "
            "reactants side must use the same compartment (unless this is an intermembrane surface reaction).";
      }
      surf_comp_id = surf_compartments[0];
    }

    bool surf_reac_has_compartment =
        surf_comp_id != COMPARTMENT_ID_INVALID && surf_comp_id != COMPARTMENT_ID_NONE;

    string comp_name;
    if (surf_reac_has_compartment) {
      comp_name = bng_data.get_compartment(surf_comp_id).name;
    }
    else {
      comp_name = compartment_id_to_str(surf_comp_id);
    }

    // do all the surface products use the same compartment?
    if (!r.is_intermembrane_surf_rxn()) {
      for (const Cplx& prod: r.products) {
        if (prod.is_surf() && prod.get_primary_compartment_id() != surf_comp_id) {

          return "Reaction rule " + r.to_str(false, false, false) + ": all compartments of surface molecules on the " +
              "products side must use the same compartment " + comp_name + " as reactants " +
              "(unless this is an intermembrane surface reaction).";
        }
      }
    }

    assert(
        (!surf_reac_has_compartment ||
         !bng_data.get_compartment(surf_comp_id).is_3d) &&
        "Should have been checked when compartments are assigned to reactants/products");
  }
  else {
    // uses_in_out_compartment
    // 3. In a surface reaction,
    //    b. all volume substances must have one of @IN and @OUT compartments
    //    and all surface substances must not have a compartment
    assert(r.is_bimol());
    assert(r.reactants[0].is_surf() == r.reactants[1].is_vol() && "There must be 1 vol and 1 surf reactant");
    CHECK(check_vol_have_in_out_compartment(r, r.reactants));
    CHECK(check_vol_have_in_out_compartment(r, r.products));
    CHECK(check_surf_have_no_compartment(r, r.reactants));
    CHECK(check_surf_have_no_compartment(r, r.products));
  }

  return "";
}


static void set_substances_orientation_from_in_out_compartments(CplxVector& substances) {
  for (Cplx& subst: substances) {
    if (subst.is_vol()) {
      if (subst.get_primary_compartment_id() == COMPARTMENT_ID_IN) {
        subst.set_orientation(ORIENTATION_DOWN);
      }
      else if (subst.get_primary_compartment_id() == COMPARTMENT_ID_OUT) {
        subst.set_orientation(ORIENTATION_UP);
      }
      else {
        assert(false);
      }
    }
  }
}


static std::string set_vol_rxn_substance_orientation_from_compartment(
    const BNGData& bng_data, RxnRule& r, const compartment_id_t surf_comp_id, Cplx& substance) {

  assert(substance.is_vol());
  assert(!is_in_out_compartment_id(surf_comp_id));

  compartment_id_t vol_comp_id = substance.get_primary_compartment_id();
  assert(!is_in_out_compartment_id(vol_comp_id) && "Handled elsewhere");

  if (surf_comp_id != COMPARTMENT_ID_NONE) {
    const Compartment& surf_comp = bng_data.get_compartment(surf_comp_id);

    if (vol_comp_id != COMPARTMENT_ID_NONE) {
      const Compartment& vol_comp = bng_data.get_compartment(vol_comp_id);

      // up or down?
      if (surf_comp.parent_compartment_id == vol_comp.id) {
        substance.set_orientation(ORIENTATION_UP);
      }
      else if (surf_comp.children_compartments.count(vol_comp.id)) {
        substance.set_orientation(ORIENTATION_DOWN);
      }
      else {
        return "Reaction rule " + r.to_str(false, true, false) + ": surface reactants use compartment " + surf_comp.name +
            " but volume reactant's compartment " + vol_comp.name + " is neither parent nor child of it.";
      }
    }
    else {
      // no volume compartment set
      substance.set_orientation(ORIENTATION_NONE); // ANY
    }
  }
  else {
    // surface compartment is not known
    if (vol_comp_id != COMPARTMENT_ID_NONE) {
      const Compartment& vol_comp = bng_data.get_compartment(vol_comp_id);
      return "Reaction rule " + r.to_str(false, true, false) + ": surface reactants do not have a compartment specified " +
          "but volume reactant's compartment " + vol_comp.name + " is set and it is not clear what whether " +
          "it shoudl be reacting from inside or outside, either specify a surface compartment for the surface reactant or use @IN or @OUT " +
          "compartment class.";
    }
    else {
      // no volume compartment set
      substance.set_orientation(ORIENTATION_NONE); // ANY
    }
  }

  return "";
}


// returns empty string if there was no error, otherwise the returned string specifies
// error message
std::string check_compartments_and_set_orientations(const BNGData& bng_data, RxnRule& r) {
  assert(r.is_finalized() && "Types of substances (vol,surf,...) are not set without finalization.");

  compartment_id_t surf_comp_id;
  bool uses_in_out_compartments;
  bool has_surf_mols_and_vol_compartments;
  CHECK(
      check_surface_compartments(
          bng_data, r, surf_comp_id, uses_in_out_compartments, has_surf_mols_and_vol_compartments)
  );

  if (surf_comp_id == COMPARTMENT_ID_NONE && !uses_in_out_compartments && !has_surf_mols_and_vol_compartments) {
    // does not have surface reactants, nothing to do,
    // this is for example the case when MCell style orientations are specified by the user
    return "";
  }

  if (uses_in_out_compartments) {
    // with IN/OUT compartments the orientation assignment is easy - IN is DOWN and OUT is UP
    set_substances_orientation_from_in_out_compartments(r.reactants);
    set_substances_orientation_from_in_out_compartments(r.products);

    // we also must set orientation of the surface reactant,
    // default orientation of a surface reactant when @IN or @OUT is set is UP
    for (Cplx& c: r.reactants) {
      if (c.is_surf()) {
        assert(c.get_orientation() == ORIENTATION_NONE || c.get_orientation() == ORIENTATION_UP);
        // this also controls the probability factor given to this reaction,
        // see RxnClass::compute_pb_factor
        c.set_orientation(ORIENTATION_UP);
      }
    }
  }
  else {
    // set orientations of reactants
    for (Cplx& reac: r.reactants) {
      if (reac.is_surf()) {
        // overwrite orientation if compartments are used or it was not set
        if (reac.has_compartment() || reac.get_orientation() == ORIENTATION_NONE) {
          reac.set_orientation(ORIENTATION_UP);
        }
      }
      else if (reac.is_vol()) {
        CHECK(set_vol_rxn_substance_orientation_from_compartment(bng_data, r, surf_comp_id, reac));
      }
      else {
        assert(reac.is_reactive_surface());
      }
    }

    // set orientations of products
    for (Cplx& prod: r.products) {
      if (prod.is_surf()) {
        if (prod.has_compartment() || prod.get_orientation() == ORIENTATION_NONE) {
          prod.set_orientation(ORIENTATION_UP);
        }
      }
      else if (prod.is_vol()) {
        CHECK(set_vol_rxn_substance_orientation_from_compartment(bng_data, r, surf_comp_id, prod));
      }
      else {
        // TODO: there should be a general check somewhere for this
        return "Reaction rule " + r.name + ": reactive surface " + prod.to_str() + " cannot be a product.";
      }
    }
  }
  return "";
}

} // namespace BNG
