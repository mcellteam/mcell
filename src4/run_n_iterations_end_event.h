/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
