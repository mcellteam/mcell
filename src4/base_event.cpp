/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <iostream>

#include "base_event.h"

using namespace std;

namespace MCell {

void BaseEvent::dump(const std::string ind) const {
  cout << ind << "event_time: \t\t" << event_time << " [double] \t\t\n";
  cout << ind << "periodicity_interval: \t\t" << periodicity_interval << " [double] \t\t\n";
  cout << ind << "type_index: \t\t" << type_index << " [event_type_index_t] \t\t\n";
}

void BaseEvent::to_data_model(Json::Value& mcell_node) const {
  // does nothing by default
}

} // namespace mcell
