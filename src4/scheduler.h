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
  // scheduler becomes owner of the base_event object
  void schedule_event(BaseEvent* event);

  // returns the time of next event
  float_t get_next_event_time();

  // returns time of the event that was handled
  EventExecutionInfo handle_next_event();

  void dump() const;

  void to_data_model(Json::Value& mcell_node) const;

private:
  Calendar calendar;
};

} // namespace mcell

#endif // SRC4_SCHEDULER_H_
