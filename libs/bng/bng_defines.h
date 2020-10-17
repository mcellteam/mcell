// this file contains BNG-specific defines

#ifndef __BNG_DEFINES_H__
#define __BNG_DEFINES_H__

#include <stdint.h>
#include <string>
#include <set>
#include <iostream>

#include "bng/defines_shared.h"

namespace BNG {

using BNGCommon::float_t;
using BNGCommon::EPS;
using BNGCommon::BNG_PI;
using BNGCommon::BNG_N_AV;
using BNGCommon::FLT_GIGANTIC;
using BNGCommon::fabs_f;
using BNGCommon::cmp_eq;
using BNGCommon::distinguishable_f;
using BNGCommon::sqrt_f;
using BNGCommon::pow_f;
using BNGCommon::floor_f;
using BNGCommon::round_f;

// arbitrary limit that is checked when we are computing all the product sets for
// a given reaction rule, only to show to the user that something is probably wrong
// (and we cannot deal with huge numbers yet, although the 1024 is still ok)
const uint MAX_PRODUCT_SETS_PER_RXN = 1024*1024;

// if the count of products after applying a single rxn is greater ahan this value
// and we know the count of the product, do not compute the products immediatelly
const uint MAX_IMMEDIATELLY_COMPUTED_PRODUCT_SETS_PER_RXN = 8;

const char* const MCELL_REDEFINE_PREFIX = "MCELL_REDEFINE_";

typedef uint state_id_t;
const state_id_t STATE_ID_INVALID = UINT32_MAX;

// for components that have a single state or
// for reactions that do not care about the state of the component
const state_id_t STATE_ID_DONT_CARE = UINT32_MAX - 2;

typedef uint component_type_id_t;
const component_type_id_t COMPONENT_TYPE_ID_INVALID = ID_INVALID;

typedef uint bond_value_t;
const bond_value_t BOND_VALUE_INVALID = UINT32_MAX; // same as ID_INVALID
const bond_value_t BOND_VALUE_UNBOUND = UINT32_MAX - 1;
const bond_value_t BOND_VALUE_BOUND = UINT32_MAX - 2; // for '+' in patterns such as a!+
const bond_value_t BOND_VALUE_ANY = UINT32_MAX - 3; // for '?' in patterns such as a!?


typedef uint mol_type_id_t;
const mol_type_id_t MOL_TYPE_ID_INVALID = ID_INVALID;

typedef uint compartment_id_t;
const compartment_id_t COMPARTMENT_ID_INVALID = ID_INVALID;
// - used for molecules, means that the molecule isn't in any of the named compartments
// - or we are not tracking compartments for this molecule because there are no reactions
//   for this molecule that use a compartment
const compartment_id_t COMPARTMENT_ID_NONE = ID_INVALID2;
// same as don't care, used in rxn classes
const compartment_id_t COMPARTMENT_ID_ANY = ID_INVALID3;

static std::string compartment_id_to_str(const compartment_id_t id) {
  switch (id) {
    case COMPARTMENT_ID_INVALID: return "INVALID";
    case COMPARTMENT_ID_NONE: return "NONE";
    case COMPARTMENT_ID_ANY: return "ANY";
    default: return std::to_string(id);
  }
}

typedef std::set<compartment_id_t> CompartmentIdSet;

// rxn rules are always global and presumed to be constant
typedef uint rxn_rule_id_t;
const rxn_rule_id_t RXN_RULE_ID_INVALID = ID_INVALID;

const char* const RXN_REPORT_PREFIX = "rxn_report_";
const char* const SPECIES_REPORT_PREFIX = "species_report_";
const char* const WARNINGS_REPORT_PREFIX = "warnings_report_";
const char* const REPORT_EXT = ".txt";
}

#endif // __BNG_DEFINES_H__
