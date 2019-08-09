
#include "defines.h"
#include "world_constants.h"

#include <iostream>
#include <sstream>

using namespace std;

namespace mcell {

void world_constants_t::dump() {
  cout << "time_unit: \t\t" << time_unit << " [float_t] \t\t\n";
  cout << "length_unit: \t\t" << length_unit << " [float_t] \t\t\n";
  cout << "rx_radius_3d: \t\t" << rx_radius_3d << " [float_t] \t\t\n";
  cout << "partition_edge_length: \t\t" << partition_edge_length << " [float_t] \t\t\n";
  cout << "subpartitions_per_partition_dimension: \t\t" << subpartitions_per_partition_dimension << " [uint] \t\t\n";
  cout << "subpartition_edge_length: \t\t" << subpartition_edge_length << " [float_t] \t\t\n";
}

}
