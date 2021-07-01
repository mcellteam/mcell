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

#ifndef SRC4_CUSTOM_FUNCTION_CALL_EVENT_H_
#define SRC4_CUSTOM_FUNCTION_CALL_EVENT_H_

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
class CustomFunctionCallEvent: public BaseEvent {
public:
  typedef void (*CalledFunctionType)(double, T);

  CustomFunctionCallEvent(
      CalledFunctionType function_ptr_,
      const T& function_arg_, // makes a shallow copy
      const event_type_index_t event_type_index_ = EVENT_TYPE_INDEX_CALL_END_ITERATION
  )
    : BaseEvent(event_type_index_),
      function_ptr(function_ptr_),
      function_arg(function_arg_),
      return_from_run_n_iterations(false){
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

#endif // SRC4_CUSTOM_FUNCTION_CALL_EVENT_H_
