/*
 * bng_engine.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_BNG_ENGINE_H_
#define LIBS_BNG_BNG_ENGINE_H_

#include "bng_defines.h"

#include "component.h"
#include "complex_species.h"

namespace BNG {

// Data shared among all instances of BNGEngines
// Usually constant, initialized when BNGL is parsed
class BNGData {

  std::vector<ComponentType> component_types;
  component_type_id_t next_component_type_id;


  // indexed with state_id_t
  std::vector<std::string> state_names;
};



// Each partition will have its own BNG engine,
// the contents might change a lot during execution and string comparison
// on using BNGL components we still have a way how to unify all the instances
class BNGEngine {
public:


  //component_index_t next_component_index;

  // checks cache - 2 types of caches can/cannot react
  // if not found, we need to decide
  void can_react(
      complex_species_index_t a,
      complex_species_index_t b,
      small_vector<complex_species_index_t>& reactions
  );

  void get_reactants(
      complex_species_index_t b,
      small_vector<complex_species_index_t>& reactions
  );

  // search whether two molecules can react is done
  std::vector<ComplexSpecies> complex_species;


  // cache of complex species indices that can interact together

  const BNGData& bng_data;
};



} /* namespace BNG */

#endif /* LIBS_BNG_BNG_ENGINE_H_ */
