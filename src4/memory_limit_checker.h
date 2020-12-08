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

#ifndef SRC4_MEMORY_LIMIT_CHECKER_H_
#define SRC4_MEMORY_LIMIT_CHECKER_H_

#include "defines.h"
#include "libs/cpptime/cpptime.h"

namespace MCell {

class World;

class MemoryLimitChecker {
public:
  MemoryLimitChecker() :
    world(nullptr),
    limit_gb(-1),
    exit_when_over_limit(false),
    over_limit(false),
    timer(nullptr), created_timer_id(0) {
  }

  ~MemoryLimitChecker();

  // does nothing when limit_in_gb_ <= 0
  void start_timed_check(
      World* world_, const int limit_gb_, const bool exit_when_over_limit_ = true);

  void stop_timed_check();

  bool is_over_memory_limit() const {
    return over_limit;
  }

private:
  // set in start_timed_check
  World* world;
  int limit_gb; // -1 to ignore
  bool exit_when_over_limit;

  bool over_limit; // safe to read asynchronously

  CppTime::Timer* timer;
  CppTime::timer_id created_timer_id;
};

} /* namespace MCell */

#endif /* SRC4_MEMORY_LIMIT_CHECKER_H_ */
