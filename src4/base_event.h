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
const event_type_index_t EVENT_TYPE_INDEX_PERIODIC_CALL = 0;
const event_type_index_t EVENT_TYPE_INDEX_RELEASE = 200;
// first counting and visualization output is done after release
const event_type_index_t EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT = 290;
// viz_output is special in the way that simulation is terminated herein the last iteration
const event_type_index_t EVENT_TYPE_INDEX_VIZ_OUTPUT = 300;

// simulation end check event is scheduled for each iteration,
// this allows us to safely terminate after all observables were collected
// (counts and viz output must be before diffuse&react for MCell3 compatibility)
const event_type_index_t EVENT_TYPE_INDEX_SIMULATION_END_CHECK = 301;

const event_type_index_t EVENT_TYPE_INDEX_RXN_CLASS_CLEANUP = 470;
const event_type_index_t EVENT_TYPE_INDEX_SPECIES_CLEANUP = 480;

const event_type_index_t EVENT_TYPE_INDEX_SORT_MOLS_BY_SUBPART = 490;
const event_type_index_t EVENT_TYPE_INDEX_DIFFUSE_REACT = 500;  // this event spans the whole time step
const event_type_index_t EVENT_TYPE_INDEX_DEFRAGMENTATION = 900;


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
  virtual float_t get_secondary_ordering_value() { return 0; }

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
  virtual float_t get_max_time_up_to_next_barrier() const {
    // the subclass must return true in may_be_blocked_by_barrier_and_needs_set_time_step
    assert(false && "Only overridden variant of this method may be called.");
    return 0;
  }

  // if an event is blocked
  virtual void set_barrier_time_for_next_execution(const float_t time_step) {
    // the subclass must return true in may_be_blocked_by_barrier_and_needs_set_time_step
    assert(false && "Only overridden variant of this method may be called.");
  }

  // time when this object's step() method will be called
  float_t event_time;

  // once this event is executed, schedule next one after this interval
  // do not schedule if the value is 0
  float_t periodicity_interval;

  // this value specifies both id of the event and ordering when multiple
  // events of a different type are sheduled for the same time
  event_type_index_t type_index;
};

}

#endif // SRC4_BASE_EVENT_H_
