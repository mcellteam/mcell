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

#ifndef SRC4_BASE_EVENT_H_
#define SRC4_BASE_EVENT_H_

#include "defines.h"

namespace Json {
class Value;
}

namespace MCell {

class World;

typedef int event_type_index_t;

// Value specifies ordering when two events are scheduled for exactly the same time
// The 'holes' are there on purpose for ordering of external events
const event_type_index_t EVENT_TYPE_INDEX_INVALID = -1;

// must be the very first event in an iteration
const event_type_index_t EVENT_TYPE_INDEX_CALL_START_ITERATION_CHECKPOINT = 0;

const event_type_index_t EVENT_TYPE_INDEX_RELEASE = 200;
// first counting and visualization output is done after release
const event_type_index_t EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT = 290;
// viz_output is special in the way that simulation is terminated herein the last iteration
const event_type_index_t EVENT_TYPE_INDEX_VIZ_OUTPUT = 300;

// simulation end check event is scheduled for each iteration,
// this allows us to safely terminate after all observables were collected
// (counts and viz output must be before diffuse&react for MCell3 compatibility)
const event_type_index_t EVENT_TYPE_INDEX_SIMULATION_END_CHECK = 301;

const event_type_index_t EVENT_TYPE_INDEX_RXN_CLASS_CLEANUP = 400;
const event_type_index_t EVENT_TYPE_INDEX_SPECIES_CLEANUP = 410;
const event_type_index_t EVENT_TYPE_INDEX_SORT_MOLS_BY_SUBPART = 420;

const event_type_index_t EVENT_TYPE_INDEX_CLAMP_RELEASE = 490;
const event_type_index_t EVENT_TYPE_INDEX_DIFFUSE_REACT = 500;  // this event spans the whole time step
const event_type_index_t EVENT_TYPE_INDEX_DEFRAGMENTATION = 900;
const event_type_index_t EVENT_TYPE_INDEX_MOL_SHUFFLE = 910;

const event_type_index_t EVENT_TYPE_INDEX_BARRIER = 980;
const event_type_index_t EVENT_TYPE_INDEX_CALL_END_ITERATION = 990;

/**
 * Base class for all events.
 * Should be independent on the mcell world in order to
 * integrate also other simulators.
 */
class BaseEvent {
public:
  BaseEvent(const event_type_index_t t) :
    event_time(TIME_INVALID), periodicity_interval(0), type_index(t) {
  }
  virtual ~BaseEvent() {};
  virtual void step() = 0;
  virtual void dump(const std::string ind = "") const;
  virtual void to_data_model(Json::Value& mcell_node) const;

  // some events such as release events have their event time set for
  // the beginning of a timestep but internally they need to be ordered
  // also according to another value such as actual release time
  virtual bool needs_secondary_ordering() { return false; }
  virtual double get_secondary_ordering_value() { return 0; }

  // if this event should be rescheduled, updates event_time and
  // possibly other attributes,
  // returns true if even should be rescheduled, false if event
  // should be removed from the schedule
  virtual bool update_event_time_for_next_scheduled_time() {
    // handling the simple case when periodicity_interval is 0 or not
    if (periodicity_interval == 0) {
      return false;
    }
    else {
      event_time = event_time + periodicity_interval;
      return true;
    }
  }

  // - events such as VizOutput, Count, or Release are barriers
  //   to events whose execution spans over a certain time step
  //   such as DiffuseReactEvent
  // - is_blocked_by_barrier_and_needs_set_time_step must return false
  //   when is_barrier returns true
  virtual bool is_barrier() const { return false; }

  // - some events represent simulation over a time step such as
  //   DiffuseReactEvent, the maximum timestep for such events
  //   must be limited so that the barrier event (such as molecule count)
  //   has data valid for its time
  // - is_barrier must return false
  //   when is_blocked_by_barrier_and_needs_set_time_step returns true
  // - periodicity_interval must not be 0 when this function return true
  virtual bool may_be_blocked_by_barrier_and_needs_set_time_step() const { return false; }

  // maximum search time for a barrier, i.e. we do not care whether
  // there is a barrier after the time interval returned by this method
  virtual double get_max_time_up_to_next_barrier() const {
    // the subclass must return true in may_be_blocked_by_barrier_and_needs_set_time_step
    assert(false && "Only overridden variant of this method may be called.");
    return 0;
  }

  // if an event is blocked
  virtual void set_barrier_time_for_next_execution(const double time_step) {
    // the subclass must return true in may_be_blocked_by_barrier_and_needs_set_time_step
    assert(false && "Only overridden variant of this method may be called.");
  }

  // used by checkpointing
  virtual bool return_from_run_n_iterations_after_execution() const {
    return false;
  }

  // time when this object's step() method will be called
  double event_time;

  // once this event is executed, schedule next one after this interval
  // do not schedule if the value is 0
  double periodicity_interval;

  // this value specifies both id of the event and ordering when multiple
  // events of a different type are sheduled for the same time
  event_type_index_t type_index;
};

}

#endif // SRC4_BASE_EVENT_H_
