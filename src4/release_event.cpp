

#include <iostream>

#include "release_event.h"

using namespace std;

namespace mcell {

void release_event_t::dump(const std::string ind) {
	cout << ind << "location: \t\t" << location << " [vec3_t] \t\t\n";
	cout << ind << "species_id: \t\t" << species_id << " [species_id_t] \t\t\n";
	cout << ind << "release_number: \t\t" << release_number << " [uint32_t] \t\t\n";
	std::cout << ind << "name: \t\t" << name << " [string] \t\t\n";
}

} // namespace mcell


