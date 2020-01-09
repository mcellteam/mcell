#ifndef __BNG_DEFINES_H__
#define __BNG_DEFINES_H__

#include <stdint.h>
#include <string>
#include "common_defines.h"
#include "../../src4/defines.h"

namespace BNG {

typedef Common::float_t float_t;

typedef uint state_id_t;
const state_id_t STATE_ID_INVALID = ID_INVALID;

// for components that have a single state or
// for reactions that do not care about the state of the component
const state_id_t STATE_ID_ANY = 0;

typedef uint component_type_id_t;
const component_type_id_t COMPONENT_TYPE_ID_INVALID = ID_INVALID;

//typedef uint component_index_t;
//const component_index_t COMPONENT_INDEX_INVALID = ID_INVALID;


typedef uint complex_species_index_t;
const complex_species_index_t COMPONENT_SPECIES_INDEX_INVALID = ID_INVALID;

typedef uint complex_instance_index_t;
const complex_instance_index_t COMPLEX_INSTANCE_INDEX_INVALID = ID_INVALID;

}

#endif // __BNG_DEFINES_H__
