/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#include <chrono>
#ifndef _WIN64
#include <sys/resource.h>
#endif

#include "memory_limit_checker.h"
#include "world.h"

#include "libs/gperftools/src/gperftools/malloc_extension.h"

using namespace CppTime;
using namespace std;
using namespace std::chrono;

namespace MCell {

const int PERIODICITY_SECONDS = 30;
const int KB_IN_GB = 1024*1024;


MemoryLimitChecker::~MemoryLimitChecker() {
  if (timer != nullptr) {
    timer->remove(created_timer_id);
    delete timer;
  }
}


void MemoryLimitChecker::start_timed_check(
    World* world_, const int limit_gb_, const bool exit_when_over_limit_) {

  if (limit_gb_ < 0) {
    return;
  }

  release_assert(timer == nullptr && "May be called only once");

  // not completely sure how lambdas work, copying arguments
  // to be sure that they are available even when we return from function
  world = world_;
  limit_gb = limit_gb_;
  exit_when_over_limit = exit_when_over_limit_;

#ifdef PROFILE_MEMORY
  static int counter = 0;
#endif

  timer = new Timer();
  created_timer_id = timer->add(
      seconds(0),
      [this](CppTime::timer_id) {

#ifdef PROFILE_MEMORY
          const std::string outfile = "mem" + std::to_string(counter) + ".dump";
          counter++;
          std::string data;
          MallocExtension::instance()->GetHeapSample(&data);
          ofstream out;
          out.open(outfile);
          out.write(data.c_str(), data.size());
          cout << "Written mem profile to " << outfile << "\n";
          out.close();
#endif

          uint64_t usage = get_mem_usage(); // returns value in kB

          if ((long long)usage > this->limit_gb * KB_IN_GB) {

            // skip if we are currently executing MolRxnCountEvent because this can mean
            // that flushed buffers are in an inconsistent state, this should be rare
            if (exit_when_over_limit &&
                (world->scheduler.get_event_being_executed() == nullptr ||
                 world->scheduler.get_event_being_executed()->type_index != EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT
                ) ) {

              errs() <<
                  "Memory limit of " << this->limit_gb << " GB reached, currently using " <<
                  usage / KB_IN_GB << " GB. Flushing count observable buffers and terminating simulation.\n";
              this->world->flush_buffers();
              exit(1);
            }
            else {
              over_limit = true;
            }

          }
      },
      seconds(PERIODICITY_SECONDS)
  );
}


void MemoryLimitChecker::stop_timed_check() {
  if (timer != nullptr) {
    timer->remove(created_timer_id);
    delete timer;
    timer = nullptr;
  }
}

} /* namespace MCell */
