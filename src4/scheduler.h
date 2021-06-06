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
#include <mutex>

#include "base_event.h"

namespace Json {
class Value;
}

namespace MCell {

// we should represent the time interval with a precisely
// represented floating point value
const double BUCKET_TIME_INTERVAL = 1;

class Bucket {
public:
  Bucket(double start_time_) :
    start_time(start_time_) {
  }
  ~Bucket();
  void insert(BaseEvent* event);

  void dump() const;

  double start_time;
  std::list<BaseEvent*> events;
};


typedef std::deque<Bucket> BucketDeque;


class Calendar {
public:
  Calendar() :
    cached_next_barrier_time(TIME_INVALID) {
    // create at least one item?
    queue.push_back(Bucket(TIME_SIMULATION_START));
  }
  ~Calendar() {
    // implicitly calls destructors of items in queue and
    // deletes all events
  }

  void insert(BaseEvent* event);
  double get_next_time();
  BaseEvent* pop_next();

  void dump() const;
  void to_data_model(Json::Value& mcell_node) const;
  const BaseEvent* find_next_event_with_type_index(
      const event_type_index_t event_type_index) const;
  void get_all_events_with_type_index(
      const event_type_index_t event_type_index, std::vector<BaseEvent*>& events);
  void get_all_events_with_type_index(
      const event_type_index_t event_type_index, std::vector<const BaseEvent*>& events) const;

  double get_time_up_to_next_barrier(const double current_time, const double max_time_step);

  void print_periodic_stats() const {
    std::cout << "Calendar: queue.size() = " << queue.size() << "\n";
  }

private:
  // differs from get_next_time - does not clear empty buckets
  // and returns bucket start time, not the next event time
  double get_first_bucket_start_time() {
    assert(queue.size() != 0);
    return queue.front().start_time;
  }

  double event_time_to_bucket_start_time(const double time) {
    // flooring to a multiple of BUCKET_TIME_INTERVAL
    return floor_to_multiple_f(time, BUCKET_TIME_INTERVAL);
  }

  BucketDeque::iterator get_or_create_bucket(const double time);

  void clear_empty_buckets();

  // queue might be empty
  BucketDeque queue;

  // used in get_time_up_to_next_barrier
  double cached_next_barrier_time;
};

// Structure used to return information about the event that was just handled
struct EventExecutionInfo {
  EventExecutionInfo(
      const double time_, const event_type_index_t type_index_, const bool return_from_run_iterations_) :
      time(time_), type_index(type_index_), return_from_run_iterations(return_from_run_iterations_) {
  }
  double time;
  event_type_index_t type_index;
  bool return_from_run_iterations;
};

class Scheduler {
public:
  Scheduler() :
    event_being_executed(nullptr),
    async_event_queue_lock(mtx, std::defer_lock),
    have_async_events_to_schedule(false) {
  }

  // - scheduler becomes owner of the base_event object
  // - event's time must be valid and not be in the past
  void schedule_event(BaseEvent* event);

  // - similar as schedule_event only guarded by a critical section and
  //   inserts the event into a separate queue of events,
  // - every method whose outcome might be changed by the events in the async queue
  //   must call schedule_events_from_async_queue to schedule these events correctly
  // - events that have event_time == TIME_INVALID will get time for the next iteration
  //   so that they are executed in the right order
  void schedule_event_asynchronously(BaseEvent* event);

  // returns the time of next event
  double get_next_event_time(const bool skip_async_events_check = false);

  // returns time of the event that was handled
  EventExecutionInfo handle_next_event();

  // skip events for checkpointing,
  // may take long time if periodic events are scheduled
  void skip_events_up_to_time(const double start_time);

  void dump() const;

  void to_data_model(Json::Value& mcell_node) const;

  const BaseEvent* find_next_event_with_type_index(
      const event_type_index_t event_type_index) {
    schedule_events_from_async_queue();

    return calendar.find_next_event_with_type_index(event_type_index);
  }

  void get_all_events_with_type_index(
      const event_type_index_t event_type_index, std::vector<BaseEvent*>& events) {
    schedule_events_from_async_queue();

    return calendar.get_all_events_with_type_index(event_type_index, events);
  }

  BaseEvent* get_event_being_executed() {
    return event_being_executed;
  }

  const BaseEvent* get_event_being_executed() const {
    return event_being_executed;
  }

  void print_periodic_stats() const {
    calendar.print_periodic_stats();
  }

private:
  // checks if there are events in the async queue and schedules them
  // not really const but maintains the consistency of events visible to the outside
  void schedule_events_from_async_queue();


  Calendar calendar;

  // callback might need to query which event is being executed right now
  // is nullptr when no event is running
  BaseEvent* event_being_executed;

  // lock for asynchronous scheduling of events
  std::mutex mtx;
  std::unique_lock<std::mutex> async_event_queue_lock;
  std::vector<BaseEvent*> async_event_queue;
  volatile bool have_async_events_to_schedule; // unguarded flag for fast check whether async_event_queue is empty
};

} // namespace mcell

#endif // SRC4_SCHEDULER_H_
