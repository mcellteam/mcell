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
 *
 * Made a template to hold its context as function_arg and
 * to be able to destroy it afterwards.
 */
template<typename T>
class EndIterationCallEvent: public BaseEvent {
public:
  typedef void (*CalledFunctionType)(float_t, T);

  EndIterationCallEvent(
      CalledFunctionType function_ptr_,
      const T& function_arg_, // makes a shallow copy
      const bool return_from_run_n_iterations_ = false
  )
    : BaseEvent(EVENT_TYPE_INDEX_END_ITERATION_CALL),
      function_ptr(function_ptr_),
      function_arg(function_arg_),
      return_from_run_n_iterations(return_from_run_n_iterations_){
  }

  // pointer to a function to be periodically called
  // first argument is event time, second is
  CalledFunctionType function_ptr;

  T function_arg;

  bool return_from_run_n_iterations;

  void step() override {
    assert(function_ptr != nullptr);
    function_ptr(event_time, function_arg);
  }

  virtual bool return_from_run_n_iterations_after_execution() const {
    return return_from_run_n_iterations;
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
