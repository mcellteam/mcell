/******************************************************************************
 *
 * Copyright (C) 2020 by
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

#ifndef SRC4_COUNT_BUFFER_H_
#define SRC4_COUNT_BUFFER_H_

#include <fstream>
#include "defines.h"

namespace MCell {

class CountItem {
public:
  CountItem()
    : time(TIME_INVALID), value(0) {
  }

  CountItem(const float_t time_, const float_t value_)
    : time(time_), value(value_) {
  }

  void inc_or_dec(const int sign, const int count = 1) {
    assert(sign == 1 || sign == -1);
    value += (sign * count);
  }

  float_t time; // time is in outside units, was already precomputed for printing
  float_t value; // e.g. count

  void write(std::ostream& out) const;
};

typedef small_vector<CountItem> CountItemVector;


// might need to be a template in the future
class CountBuffer {
public:

  CountBuffer(
      const std::string filename_,
      const size_t buffer_size_,
      const bool open_for_append_) :
      filename(filename_),
      buffer_size(buffer_size_),
      open_for_append(open_for_append_) {
  }

  void add(const CountItem& item) {
    if (data.size() >= buffer_size) {
      flush();
    }
    data.push_back(item);
  }

  // close, create an empty file
  void flush_and_close();

  const std::string& get_filename() const {
    return filename;
  }

private:
  // flush buffer, open output file if needed, keep file open afterwards
  void flush();

  // open file, return false if file could not be opened and error_is_fatal is false
  bool open(bool error_is_fatal = true);

private:
  // name of the output file with the full path
  std::string filename;

  // number of items to be stored, automatically flushes afterwards
  size_t buffer_size;

  // output stream
  std::ofstream fout;

  // buffer
  CountItemVector data;

  bool open_for_append;
};

typedef std::vector<CountBuffer> CountBufferVector;

} /* namespace MCell */

#endif /* SRC4_COUNT_BUFFER_H_ */
