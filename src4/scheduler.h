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

#include "defines.h"
#include "base_event.h"

namespace mcell {

// preferrably, we should represent the time interval precisely
const float_t SCHEDULER_BUCKET_TIME_INTERVAL = 1;

template<class BUCKET_ITEM_T> // used only for diffuse_or_unimol_react_action_t
class fifo_bucket_t {
public:
  fifo_bucket_t(float_t start_time_) :
    start_time(start_time_) {
  }
  ~fifo_bucket_t() {
    // inserted items are not pointers
  }
  
  
  void insert(const BUCKET_ITEM_T& event) {
    events.push_back(event);
  }

  float_t start_time;
  std::vector<BUCKET_ITEM_T> events;
};


template<class BUCKET_ITEM_T> // used only with base_event_t*, must be a pointer
class sorted_bucket_t {
  typedef BUCKET_ITEM_T bucket_event_t;
public:
  sorted_bucket_t(float_t start_time_) :
    start_time(start_time_) {
  }
  ~sorted_bucket_t() {
    for (auto it = events.begin(); it != events.end(); it++) {
      // delete remaining events, usually there should be none
      delete *it;
    }
  }
  

  void insert(BUCKET_ITEM_T event) {
    // check right away if the event belongs to the end
    if (events.empty() || cmp_lt(events.back()->event_time, event->event_time, SCHEDULER_COMPARISON_EPS)) {
      events.push_back(event);
    }
    else {
      // go through the list and put our item to the right time and right order
      auto it = events.begin();

      // find the right time
      while (it != events.end() && cmp_lt((*it)->event_time, event->event_time, SCHEDULER_COMPARISON_EPS)) {
        it++;
      }
      // find the right ordering among events with the same event_time
      while (it != events.end() && cmp_eq((*it)->event_time, event->event_time, SCHEDULER_COMPARISON_EPS)
          && (*it)->type_index <= event->type_index) {
        it++;
      }

      events.insert(it, event);
    }
  }

  float_t start_time;
  std::list<BUCKET_ITEM_T> events;
};


template<class BUCKET_T, typename BUCKET_ITEM_T>
class calendar_t {
public:
  typedef std::deque<BUCKET_T> bucket_deque_t;

  calendar_t(const float_t bucket_time_interval_)
    : bucket_time_interval(bucket_time_interval_) {
    bucket_time_interval_rcp = 1/bucket_time_interval;
    // create at least one item
    queue.push_back(BUCKET_T(TIME_SIMULATION_START));
  }
  ~calendar_t() {
    // implicitly calls destructors of items in queue and
    // deletes all events
  }


  void insert(BUCKET_ITEM_T event, float_t event_time) {
    // event_time == base_event::event_time == diffuse_or_unimol_react_action_t::scheduled_time
    float_t bucket_start_time = event_time_to_bucket_start_time(event_time);
    if (queue.empty()) {
      // no items yet - simply create new bucket and insert our event there
      queue.push_back(BUCKET_T(bucket_start_time));
      queue.front().insert(event);
    }
    else {
      // we first need to find out whether we already have a bucket for this event
      float_t first_start_time = get_first_bucket_start_time();
      assert(bucket_start_time - first_start_time >= 0 && "cannot schedule to the past"); // some eps?
      size_t buckets_from_first = (bucket_start_time - first_start_time) * bucket_time_interval_rcp;

      if (buckets_from_first < queue.size()) {
        // bucket exists
        queue[buckets_from_first].insert(event);
      }
      else {
        // we need to create new buckets
        size_t missing_buckets = buckets_from_first - queue.size() + 1;
        float_t next_time = queue.back().start_time + bucket_time_interval;
        for (size_t i = 0; i < missing_buckets; i++) {
          queue.push_back(BUCKET_T(next_time));
          next_time += bucket_time_interval;
        }
        assert(buckets_from_first < queue.size());
        queue[buckets_from_first].insert(event);
      }
    }
  }


  // returns BUCKET_INDEX_INVALID if bucket does not exist
  uint64_t get_bucket_index_for_time(const float_t time) {
    float_t bucket_start_time = event_time_to_bucket_start_time(time);
    float_t first_start_time = get_first_bucket_start_time();
    int64_t buckets_from_first = (bucket_start_time - first_start_time) * bucket_time_interval_rcp;
    if (buckets_from_first >= 0 && buckets_from_first < (int64_t)queue.size()) {
      return buckets_from_first;
    }
    else {
      return BUCKET_INDEX_INVALID;
    }
  }


  BUCKET_T& get_bucket_with_index(const size_t index) {
    assert(index < queue.size() && "Bucket with requested index does not exist");
    return queue[index];
  }


  void pop_bucket() {
    assert(!queue.empty());
    queue.pop_front();
  }


  BUCKET_ITEM_T pop_next() {
    while (queue.front().events.empty()) {
      queue.pop_front();
    }
    BUCKET_ITEM_T next_event = queue.front().events.front();
    queue.front().events.pop_front();
    return next_event;
  }

private:
  float_t get_first_bucket_start_time() {
    assert(queue.size() != 0);
    return queue.front().start_time;
  }


  float_t event_time_to_bucket_start_time(const float_t time) {
    // flooring to a multiple of BUCKET_TIME_INTERVAL
    return floor_to_multiple(time, bucket_time_interval);
  }

  // queue might be empty
  bucket_deque_t queue;
  float_t bucket_time_interval;
  float_t bucket_time_interval_rcp;
};


class scheduler_t {
public:
  scheduler_t()
    : calendar(SCHEDULER_BUCKET_TIME_INTERVAL) {
  }

  // scheduler becomes owner of the base_event object
  void schedule_event(base_event_t* event);

  // returns time of the event that was handled
  float_t handle_next_event(bool &end_simulation);

private:
  calendar_t<sorted_bucket_t<base_event_t*>, base_event_t*> calendar;
};

} // namespace mcell

#endif // SRC4_SCHEDULER_H_
