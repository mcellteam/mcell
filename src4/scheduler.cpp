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

namespace mcell {


void scheduler_t::schedule_event(base_event_t* event) {
  calendar.insert(event, event->event_time);
}


// pop next scheduled event and run its step method
float_t scheduler_t::handle_next_event(bool &end_simulation) {

  base_event_t* event = calendar.pop_next();
  assert(event != NULL && "Empty event queue - at least end simulation event should be present");
  float_t event_time = event->event_time;

#ifdef DEBUG_SCHEDULER
  event->dump("");
#endif
  event->step();

  // schedule itself for the next period or just delete
  if (event->periodicity_interval != 0) {
    event->event_time += event->periodicity_interval;
    calendar.insert(event, event->event_time);
  }
  else {
    delete event;
  }

  end_simulation = event->type_index == EVENT_TYPE_INDEX_END_SIMULATION;
  return event_time;
}

} // namespace mcell
