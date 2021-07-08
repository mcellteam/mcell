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

#include "defines.h"

#include <iostream>
#include <sstream>

#ifndef _MSC_VER
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifdef DWITHGPERFTOOLS
// using longer path to avoid collisions
#include "install_gperftools/include/profiler.h"
#endif

using namespace std;

namespace MCell {

string Vec3::to_string() const {
  stringstream ss;
  ss << *this;
  return ss.str();
}


void Vec3::dump(const std::string extra_comment, const std::string ind) const {
  cout << ind << extra_comment << *this << "\n";
}


string Vec2::to_string() const {
  stringstream ss;
  ss << *this;
  return ss.str();
}


void Vec2::dump(const std::string extra_comment, const std::string ind) const {
  cout << ind << extra_comment << *this << "\n";
}


uint64_t get_mem_usage() {
#ifdef _WIN64
  return 0;
#else
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;

  ret = getrusage(who,&usage);

  if (ret == 0) {
    return usage.ru_maxrss;
  }
  else {
    // ignoring fail
    return 0;
  }
#endif
}


} // namespace mcell
