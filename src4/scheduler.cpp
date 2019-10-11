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

#include "scheduler.h"

namespace MCell {

void Bucket::insert(BaseEvent* event) {
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


Bucket::~Bucket() {
  for (auto it = events.begin(); it != events.end(); it++) {
    // delete remaining events, usually there should be none
    delete *it;
  }
}


// insert a new item with time event->event_time, create bucket if needed
void Calendar::insert(BaseEvent* event) {
  float_t bucket_start_time = event_time_to_bucket_start_time(event->event_time);
  if (queue.empty()) {
    // no items yet - simply create new bucket and insert our event there
    queue.push_back( Bucket(bucket_start_time) );
    queue.front().insert(event);
  }
  else {
    // we first need to find out whether we already have a bucket for this event
    float_t first_start_time = get_first_bucket_start_time();
    assert(bucket_start_time - first_start_time >= 0 && "cannot schedule to the past"); // some eps?
    size_t buckets_from_first = (bucket_start_time - first_start_time) / BUCKET_TIME_INTERVAL;

    if (buckets_from_first < queue.size()) {
      // bucket exists
      queue[buckets_from_first].insert(event);
    }
    else {
      // we need to create new buckets
      size_t missing_buckets = buckets_from_first - queue.size() + 1;
      float_t next_time = queue.back().start_time + BUCKET_TIME_INTERVAL;
      for (size_t i = 0; i < missing_buckets; i++) {
        queue.push_back(Bucket(next_time));
        next_time += BUCKET_TIME_INTERVAL;
      }
      assert(buckets_from_first < queue.size());
      queue[buckets_from_first].insert(event);
    }
  }
}


void Calendar::clear_empty_buckets() {
  while (queue.front().events.empty()) {
    queue.pop_front();
  }
}

BaseEvent* Calendar::pop_next() {
  clear_empty_buckets();
  BaseEvent* next_event = queue.front().events.front();
  queue.front().events.pop_front();
  return next_event;
}

float_t Calendar::get_next_time() {
  clear_empty_buckets();
  BaseEvent* next_event = queue.front().events.front();
  return next_event->event_time;
}

void Scheduler::schedule_event(BaseEvent* event) {
  calendar.insert(event);
}

float_t Scheduler::get_next_event_time() {
  return calendar.get_next_time();
}

// pop next scheduled event and run its step method
// current_iteration is the world object's iteration
float_t Scheduler::handle_next_event() {
  BaseEvent* event = calendar.pop_next();
  assert(event != NULL && "Empty event queue - at least end simulation event should be present");
  float_t event_time = event->event_time;

#ifdef DEBUG_SCHEDULER
  event->dump("");
#endif
  event->step();

  // schedule itself for the next period or just delete
  if (event->periodicity_interval != 0) {
    event->event_time += event->periodicity_interval;
    calendar.insert(event);
  }
  else {
    delete event;
  }

  return event_time;
}

} // namespace mcell