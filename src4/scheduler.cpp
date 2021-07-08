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

      // if we already found events with our type, check secondary ordering
      if ((*it)->type_index == event->type_index && event->needs_secondary_ordering()) {
        assert((*it)->needs_secondary_ordering());

        // do we belong in front of the current event?
        if (event->get_secondary_ordering_value() < (*it)->get_secondary_ordering_value()) {
          // yes, terminate search
          break;
        }
      }

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


void Bucket::dump() const {
  for (const BaseEvent* event: events) {
    event->dump();
  }
}


// insert a new item with time event->event_time, create bucket if needed
void Calendar::insert(BaseEvent* event) {

  // update barrier cache if needed
  if (event->is_barrier() && event->event_time < cached_next_barrier_time) {
    cached_next_barrier_time = event->event_time;
  }

  // align time to multiple of one if the value is close to it
  // required for example for custom time step where 0.1 * 10 must be equal to 1
  double rounded_time = round_f(event->event_time);
  if (cmp_eq(round_f(event->event_time), event->event_time, SCHEDULER_COMPARISON_EPS)) {
    event->event_time = rounded_time;
  }

  double bucket_start_time = event_time_to_bucket_start_time(event->event_time);
  if (queue.empty()) {
    // no items yet - simply create new bucket and insert our event there
    queue.push_back( Bucket(bucket_start_time) );
    queue.front().insert(event);
  }
  else {
    // we first need to find out whether we already have a bucket for this event
    double first_start_time = get_first_bucket_start_time();
    release_assert(bucket_start_time - first_start_time >= 0 && "cannot schedule to the past"); // some eps?
    size_t buckets_from_first = (bucket_start_time - first_start_time) / BUCKET_TIME_INTERVAL;

    if (buckets_from_first < queue.size()) {
      // bucket exists
      queue[buckets_from_first].insert(event);
    }
    else {
      // we need to create new buckets
      size_t missing_buckets = buckets_from_first - queue.size() + 1;
      double next_time = queue.back().start_time + BUCKET_TIME_INTERVAL;
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


double Calendar::get_next_time() {
  clear_empty_buckets();
  assert(!queue.empty() && !queue.front().events.empty());
  BaseEvent* next_event = queue.front().events.front();
  return next_event->event_time;
}


void Calendar::dump() const {
  for (const Bucket& bucket: queue) {
    bucket.dump();
  }
}


void Calendar::to_data_model(Json::Value& mcell_node) const {
  for (const Bucket& bucket: queue) {
    for (const BaseEvent* event: bucket.events) {
      event->to_data_model(mcell_node);
    }
  }
}


const BaseEvent* Calendar::find_next_event_with_type_index(
    const event_type_index_t event_type_index) const {

  for (const Bucket& bucket: queue) {
    for (const BaseEvent* event: bucket.events) {
      if (event->type_index == event_type_index) {
        return event;
      }
    }
  }
  return nullptr;
}


void Calendar::get_all_events_with_type_index(
    const event_type_index_t event_type_index,
    std::vector<BaseEvent*>& events
) {

  for (const Bucket& bucket: queue) {
    for (BaseEvent* event: bucket.events) {
      if (event->type_index == event_type_index) {
        events.push_back(event);
      }
    }
  }
}


void Calendar::get_all_events_with_type_index(
    const event_type_index_t event_type_index,
    std::vector<const BaseEvent*>& events
) const {

  for (const Bucket& bucket: queue) {
    for (BaseEvent* event: bucket.events) {
      if (event->type_index == event_type_index) {
        events.push_back(event);
      }
    }
  }
}


// returns max_time_step if no barrier is scheduled for interval
// current_time .. current_time+max_time_step
// if such a barrier exists, returns barrier time - current_time
double Calendar::get_time_up_to_next_barrier(
    const double current_time, const double max_time_step) {

  if (cached_next_barrier_time != TIME_INVALID &&
      current_time < cached_next_barrier_time
  ) {
    return std::min(cached_next_barrier_time - current_time, max_time_step);
  }

  // go through the whole schedule and find the first barrier
  cached_next_barrier_time = TIME_FOREVER;

  // expecting that there are only events that are scheduled
  // for the future
  bool found = false;
  for (const Bucket& bucket: queue) {
    for (const BaseEvent* event: bucket.events) {
      if (event->is_barrier()) {
        assert(event->event_time >= current_time);
        cached_next_barrier_time = event->event_time;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }

  return std::min(cached_next_barrier_time - current_time, max_time_step);
}


void Scheduler::schedule_event(BaseEvent* event) {
  release_assert(event->event_time != TIME_INVALID);
  calendar.insert(event);
}


void Scheduler::schedule_event_asynchronously(BaseEvent* event) {
  // may be called multiple times at the same moment e.g. when
  // adding a checkpointing event based on some timer in Python
  async_event_queue_lock.lock();
  have_async_events_to_schedule = true;
  async_event_queue.push_back(event);
  async_event_queue_lock.unlock();
}


void Scheduler::schedule_events_from_async_queue() {
  if (!have_async_events_to_schedule) {
    return;
  }

  async_event_queue_lock.lock();
  for (BaseEvent* e: async_event_queue) {

    if (e->event_time == TIME_INVALID) {
      // real asynchronous event,
      // schedule correctly for an upcoming iteration while making sure
      // that all the events from the current iteration are finished so that
      // the event_type_index is correctly followed
      e->event_time = floor_f(get_next_event_time(true) + 1);
    }
    schedule_event(e);
  }
  async_event_queue.clear();
  have_async_events_to_schedule = false;
  async_event_queue_lock.unlock();
}


double Scheduler::get_next_event_time(const bool skip_async_events_check) {
  if (!skip_async_events_check) {
    schedule_events_from_async_queue();
  }

  return calendar.get_next_time();
}


// pop next scheduled event and run its step method
EventExecutionInfo Scheduler::handle_next_event() {

  // first check if there are any
  schedule_events_from_async_queue();

  BaseEvent* event = calendar.pop_next();
  assert(event != NULL && "Empty event queue - at least end simulation event should be present");
  double event_time = event->event_time;

  if (event->may_be_blocked_by_barrier_and_needs_set_time_step()) {
    double max_time_step = calendar.get_time_up_to_next_barrier(
        event->event_time, event->get_max_time_up_to_next_barrier());
    event->set_barrier_time_for_next_execution(max_time_step);
  }

#ifdef DEBUG_SCHEDULER
  event->dump("");
#endif
  event_being_executed = event;
  event->step();
  event_being_executed = nullptr;

  event_type_index_t type_index = event->type_index;
  bool return_from_run_iterations = event->return_from_run_n_iterations_after_execution();

  // schedule itself for the next period or just delete
  double next_time;
  bool to_schedule = event->update_event_time_for_next_scheduled_time();
  if (to_schedule) {
    calendar.insert(event);
  }
  else {
    delete event;
  }

  return EventExecutionInfo(event_time, type_index, return_from_run_iterations);
}


void Scheduler::skip_events_up_to_time(const double start_time) {

  // need to deal with imprecisions, e.g. 0.0000005 * 10^6 ~= 5.0000000000008
  while (calendar.get_next_time() < start_time - EPS) {
    BaseEvent* event = calendar.pop_next();
    bool to_schedule = event->update_event_time_for_next_scheduled_time();
    if (to_schedule) {
      calendar.insert(event);
    }
    else {
      delete event;
    }
  }
}


void Scheduler::dump() const {
  calendar.dump();
}


void Scheduler::to_data_model(Json::Value& mcell_node) const {
  // go through all events and run their conversion,
  // for many events, the conversion does nothing
  calendar.to_data_model(mcell_node);
}


} // namespace mcell
