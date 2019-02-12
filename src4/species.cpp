/*
 * species.cpp
 *
 *  Created on: Feb 10, 2019
 *      Author: adam
 */
#include <iostream>

#include "species.h"

using namespace std;

namespace mcell {

void species_t::dump(const std::string ind) {
  cout << ind <<"species_id: \t\t" << species_id << " [uint16_t] \t\t/* Unique ID for this species */\n";
  cout << ind <<"mcell_species_id: \t\t" << mcell3_species_id << " [uint32_t] \t\t/* Unique ID for this species from mcell3 representation*/\n";
  cout << ind <<"name: *\t\t" << name << " [string] \t\t/* Symbol table entry (name) */\n";
  cout << ind <<"D: \t\t" << D << " [float_t] \t\t/* Diffusion constant */\n";
  cout << ind <<"space_step: \t\t" << space_step << " [float_t] \t\t/* Characteristic step length */\n";
  cout << ind <<"time_step: \t\t" << time_step << " [float_t] \t\t/* Minimum (maximum?) sensible timestep */\n";
}


} /* namespace mcell */
