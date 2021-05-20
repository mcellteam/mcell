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
