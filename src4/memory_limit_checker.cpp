/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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

              this->world->fatal_error(
                  "Memory limit of " + to_string(this->limit_gb) + " GB reached, currently using " +
                  to_string(usage / KB_IN_GB) + " GB. Flushing count observable buffers and terminating simulation.\n"
              );
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
