// this file contains BNG-specific defines

#ifndef __BNG_DEFINES_H__
#define __BNG_DEFINES_H__

#include <stdint.h>
#include <string>
#include <set>
#include <iostream>

#include "bng/shared_defines.h"

#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount __popcnt
#endif

namespace BNG {

using BNGCommon::f_to_str;
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

// if unimol rxn probability is very high, this causes the simulation to practically fails
// because the time interval between molecule creation and such unimol reaction is close to 0
const float_t MAX_UNIMOL_RXN_PROBABILITY = 1e8;

// if the count of products after applying a single rxn is greater ahan this value
// and we know the count of the product, do not compute the products immediatelly
const uint MAX_IMMEDIATELLY_COMPUTED_PRODUCT_SETS_PER_RXN = 8;

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


typedef uint elem_mol_type_id_t;
const elem_mol_type_id_t MOL_TYPE_ID_INVALID = ID_INVALID;

typedef uint compartment_id_t;
const compartment_id_t COMPARTMENT_ID_INVALID = UINT32_MAX; // same as ID_INVALID;
// - used for molecules, means that the molecule isn't in any of the named compartments
// - or we are not tracking compartments for this molecule because there are no reactions
//   for this molecule that use a compartment
// - used also as ANY compartment
const compartment_id_t COMPARTMENT_ID_NONE = UINT32_MAX - 1;
// compartment classes
const compartment_id_t COMPARTMENT_ID_IN = UINT32_MAX - 3;
const compartment_id_t COMPARTMENT_ID_OUT = UINT32_MAX - 4;

const char* const COMPARTMENT_NAME_IN = "IN";
const char* const COMPARTMENT_NAME_OUT = "OUT";

static bool is_in_out_compartment_id(const compartment_id_t id) {
  assert(id != COMPARTMENT_ID_INVALID);
  return id == COMPARTMENT_ID_IN || id == COMPARTMENT_ID_OUT;
}

static bool is_specific_compartment_id(const compartment_id_t id) {
  assert(id != COMPARTMENT_ID_INVALID);
  return id != COMPARTMENT_ID_NONE && !is_in_out_compartment_id(id);
}

// returns COMPARTMENT_ID_INVALID if name is not IN or OUT
static compartment_id_t get_in_or_out_compartment_id(const std::string& name) {
  if (name == COMPARTMENT_NAME_IN) {
    return COMPARTMENT_ID_IN;
  }
  else if (name == COMPARTMENT_NAME_OUT) {
    return COMPARTMENT_ID_OUT;
  }
  else {
    return COMPARTMENT_ID_INVALID;
  }
}

static std::string compartment_id_to_str(const compartment_id_t id) {
  switch (id) {
    case COMPARTMENT_ID_INVALID: return "INVALID";
    case COMPARTMENT_ID_NONE: return "NONE";
    case COMPARTMENT_ID_IN: return COMPARTMENT_NAME_IN;
    case COMPARTMENT_ID_OUT: return COMPARTMENT_NAME_OUT;
    default: return std::to_string(id);
  }
}

// same string are defined as API::ALL... in libmcell's gen_constants.h but
// we want BNG lib to be independent
const char* const ALL_MOLECULES = "ALL_MOLECULES";
const char* const ALL_VOLUME_MOLECULES = "ALL_VOLUME_MOLECULES";
const char* const ALL_SURFACE_MOLECULES = "ALL_SURFACE_MOLECULES";
const int NUM_GENERAL_SPECIES = 3;

static bool is_species_superclass(const std::string& name) {
  return name == ALL_MOLECULES || name == ALL_VOLUME_MOLECULES || name == ALL_SURFACE_MOLECULES;
}


typedef std::set<compartment_id_t> CompartmentIdSet;

// rxn rules are always global and presumed to be constant
typedef uint rxn_rule_id_t;
const rxn_rule_id_t RXN_RULE_ID_INVALID = ID_INVALID;

typedef uint reactant_class_id_t;
const reactant_class_id_t REACTANT_CLASS_ID_INVALID = ID_INVALID;

typedef uint_set<reactant_class_id_t> ReactantClassIdSet;

const char PATH_SEPARATOR =
#ifdef _WIN64
                            '\\';
#else
                            '/';
#endif

const char* const REPORT_DIR = "reports";
const char* const RXN_REPORT_PREFIX = "rxn_report_";
const char* const SPECIES_REPORT_PREFIX = "species_report_";
const char* const WARNINGS_REPORT_PREFIX = "warnings_report_";
const char* const REPORT_EXT = ".txt";
}

#endif // __BNG_DEFINES_H__
