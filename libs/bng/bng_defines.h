// this file contains BNG-specific defines

#ifndef __BNG_DEFINES_H__
#define __BNG_DEFINES_H__

#include <stdint.h>
#include <string>
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

typedef uint state_id_t;
const state_id_t STATE_ID_INVALID = UINT32_MAX;

// for components that have a single state or
// for reactions that do not care about the state of the component
const state_id_t STATE_ID_DONT_CARE = UINT32_MAX - 2;

typedef uint component_type_id_t;
const component_type_id_t COMPONENT_TYPE_ID_INVALID = ID_INVALID;

typedef uint bond_value_t;
const bond_value_t BOND_VALUE_INVALID = UINT32_MAX; // same as ID_INVALID
const bond_value_t BOND_VALUE_ANY = UINT32_MAX - 2; // for '+' in patterns such as a!+
const bond_value_t BOND_VALUE_NO_BOND = UINT32_MAX - 3;

typedef uint mol_type_id_t;
const mol_type_id_t MOL_TYPE_ID_INVALID = ID_INVALID;

// rxn rules are always global and presumed to be constant
typedef uint rxn_rule_id_t;
const rxn_rule_id_t RXN_RULE_ID_INVALID = ID_INVALID;

}

#endif // __BNG_DEFINES_H__
