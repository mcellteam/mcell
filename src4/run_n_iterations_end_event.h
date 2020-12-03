/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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


#ifndef SRC4_RUN_N_ITERATIONS_END_EVENT_H_
#define SRC4_RUN_N_ITERATIONS_END_EVENT_H_

#include <iostream>
#include <string>

#include "base_event.h"

namespace MCell {

/**
 * This is a dummy event that only serves as a marker
 * of time where we must check whether we already reached target number
 * of iterations.
 */
class RunNIterationsEndEvent: public BaseEvent {
public:
  RunNIterationsEndEvent()
    : BaseEvent(EVENT_TYPE_INDEX_SIMULATION_END_CHECK) {
  }

  void step() override {
    // empty
  }

  bool is_barrier() const override {
    return true;
  }

  void dump(const std::string ind) const override {
    std::cout << ind << "Simulation end check event\n";
    std::string ind2 = ind + "  ";
    BaseEvent::dump(ind2);
  }
};

} /* namespace MCell */

#endif /* SRC4_RUN_N_ITERATIONS_END_EVENT_H_ */
