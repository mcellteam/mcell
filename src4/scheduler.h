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

#ifndef SRC4_SCHEDULER_H_
#define SRC4_SCHEDULER_H_

#include <deque>
#include <list>

#include "base_event.h"

namespace Json {
class Value;
}

namespace MCell {

// we should represent the time interval with a precisely
// represented floating point value
const float_t BUCKET_TIME_INTERVAL = 1;

class Bucket {
public:
  Bucket(float_t start_time_) :
    start_time(start_time_) {
  }
  ~Bucket();
  void insert(BaseEvent* event);

  void dump() const;

  float_t start_time;
  std::list<BaseEvent*> events;
};


typedef std::deque<Bucket> BucketDeque;


class Calendar {
public:
  Calendar() {
    // create at least one item?
    queue.push_back(Bucket(TIME_SIMULATION_START));
  }
  ~Calendar() {
    // implicitly calls destructors of items in queue and
    // deletes all events
  }

  void insert(BaseEvent* event);
  float_t get_next_time();
  BaseEvent* pop_next();

  void dump() const;
  void to_data_model(Json::Value& mcell_node) const;
  const BaseEvent* find_next_event_with_type_index(
      const event_type_index_t event_type_index) const;
  void get_all_events_with_type_index(
      const event_type_index_t event_type_index, std::vector<BaseEvent*>& events);
  void get_all_events_with_type_index(
      const event_type_index_t event_type_index, std::vector<const BaseEvent*>& events) const;

  float_t get_time_up_to_next_barrier(const float_t current_time, const float_t max_time_step) const;

  void print_periodic_stats() const {
    std::cout << "Calendar: queue.size() = " << queue.size() << "\n";
  }

private:
  float_t get_first_bucket_start_time() {
    assert(queue.size() != 0);
    return queue.front().start_time;
  }

  float_t event_time_to_bucket_start_time(const float_t time) {
    // flooring to a multiple of BUCKET_TIME_INTERVAL
    return floor_to_multiple(time, BUCKET_TIME_INTERVAL);
  }

  BucketDeque::iterator get_or_create_bucket(const float_t time);

  void clear_empty_buckets();

  // queue might be empty
  BucketDeque queue;
};

// Structure used to return information about the event that was just handled
struct EventExecutionInfo {
  EventExecutionInfo(float_t time_, event_type_index_t type_index_)
    : time(time_), type_index(type_index_) {
  }
  float_t time;
  event_type_index_t type_index;
};

class Scheduler {
public:
  Scheduler() :
    event_being_executed(nullptr) {
  }

  // scheduler becomes owner of the base_event object
  void schedule_event(BaseEvent* event);

  // returns the time of next event
  float_t get_next_event_time();

  // returns time of the event that was handled
  EventExecutionInfo handle_next_event();

  // skip events for checkpointing,
  // may take long time if periodic events are scheduled
  void skip_events_up_to_time(const float_t start_time);

  void dump() const;

  void to_data_model(Json::Value& mcell_node) const;

  const BaseEvent* find_next_event_with_type_index(
      const event_type_index_t event_type_index) const {
    return calendar.find_next_event_with_type_index(event_type_index);
  }

  void get_all_events_with_type_index(
      const event_type_index_t event_type_index, std::vector<BaseEvent*>& events) {
    return calendar.get_all_events_with_type_index(event_type_index, events);
  }

  void get_all_events_with_type_index(
      const event_type_index_t event_type_index, std::vector<const BaseEvent*>& events) const {
    return calendar.get_all_events_with_type_index(event_type_index, events);
  }

  BaseEvent* get_event_being_executed() {
    return event_being_executed;
  }

  void print_periodic_stats() const {
    calendar.print_periodic_stats();
  }

private:
  Calendar calendar;

  // callback might need to query which event is being executed right now
  // is nullptr when no event is running
  BaseEvent* event_being_executed;
};

} // namespace mcell

#endif // SRC4_SCHEDULER_H_
