/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#ifndef SRC4_END_SIMULATION_EVENT_H_
#define SRC4_END_SIMULATION_EVENT_H_

#include "base_event.h"

namespace mcell {

class end_simulation_event_t: public mcell::base_event_t {
public:
  end_simulation_event_t() :
    base_event_t(EVENT_TYPE_INDEX_END_SIMULATION) {
  }
  void step() {
    // does nothing, this type of event is detected in the scheduler
  }
  void dump(const std::string indent) {
    // TODO
  }
};

} // namespace mcell

#endif // SRC4_END_SIMULATION_EVENT_H_
