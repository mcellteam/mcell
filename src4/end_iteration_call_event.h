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

#ifndef SRC4_END_ITERATION_CALL_EVENT_H_
#define SRC4_END_ITERATION_CALL_EVENT_H_

#include "base_event.h"

namespace MCell {

/**
 * General event that allows to call a function periodically.
 * Used for example to check that the user pressed ctrl-c
 * when running inside the Python interpreter.
 * Always executed at the end of an iteration.
 */
class EndIterationCallEvent: public BaseEvent {
public:
  EndIterationCallEvent(World* /*world_*/)
    : BaseEvent(EVENT_TYPE_INDEX_END_ITERATION_CALL),
      function_ptr(nullptr),
      function_arg(nullptr) {
  }

  // pointer to a function to be periodically called
  void (*function_ptr)(float_t, void*);

  void* function_arg;

  void step() override {
    assert(function_ptr != nullptr);
    function_ptr(event_time, function_arg);
  }

  void dump(const std::string ind) const override {
    std::cout << ind << "Periodic Call Event:\n";
    std::string ind2 = ind + "  ";
    BaseEvent::dump(ind2);
    std::cout << ind2 << "function_ptr:\t\t" << std::hex << (void*)function_ptr << std::dec << "\n";
  }
};

} // namespace mcell

#endif // SRC4_END_ITERATION_CALL_EVENT_H_
