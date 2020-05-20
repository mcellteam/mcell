
#include "defines_shared.h" // includes also bng_config.h

#include <iostream>

using namespace std;

namespace BNG {

void BNGNotifications::dump() const {
  cout << "BNGNotifications:\n";
  cout << "  rxn_probability_changed: \t\t" << rxn_probability_changed << " [bool] \t\t\n";
}


void BNGConfig::dump() const {
  cout << "BNGConfig:\n";
  cout << "  time_unit: \t\t" << time_unit << " [float_t] \t\t\n";
  cout << "  length_unit: \t\t" << length_unit << " [float_t] \t\t\n";
  cout << "  grid_density: \t\t" << grid_density << " [float_t] \t\t\n";
  cout << "  rx_radius_3d: \t\t" << rx_radius_3d << " [float_t] \t\t\n";
  cout << "  debug_reactions: \t\t" << debug_reactions << " [bool] \t\t\n";

  notifications.dump();
}

} // namespace BNG
