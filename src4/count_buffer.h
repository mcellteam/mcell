/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_COUNT_BUFFER_H_
#define SRC4_COUNT_BUFFER_H_

#include <fstream>
#include "defines.h"
#include "generated/gen_constants.h"

namespace MCell {

using API::CountOutputFormat;

class CountItem {
public:
  CountItem()
    : column_index(UINT_INVALID), time(TIME_INVALID), value(0) {
  }

  CountItem(const double time_, const double value_)
    : column_index(UINT_INVALID), time(time_), value(value_) {
  }

  void inc_or_dec(const int sign, const int count = 1) {
    assert(sign == 1 || sign == -1);
    value += (sign * count);
  }

  uint column_index; // index in CountBuffer
  double time; // time is in outside units, was already precomputed for printing
  double value; // e.g. count

  void write_as_dat(std::ostream& out) const;
};

typedef small_vector<CountItem> CountItemVector;


// might need to be a template in the future
class CountBuffer {
public:

  CountBuffer(
      const CountOutputFormat output_format_,
      const std::string filename_,
      const std::vector<std::string> column_names_,
      const size_t buffer_size_,
      const bool open_for_append_) :
      output_format(output_format_),
      filename(filename_),
      column_names(column_names_),
      buffer_size(buffer_size_),
      open_for_append(open_for_append_) {
    assert(output_format != CountOutputFormat::UNSET);
    assert(!column_names.empty());
    // there is a single column
    data.resize(column_names.size());
  }

  void add(const CountItem& item) {
    assert(item.column_index < data.size());
    if (data[item.column_index].size() >= buffer_size) {
      flush();
    }
    data[item.column_index].push_back(item);
  }

  // close, create an empty file
  void flush_and_close();

  const std::string& get_filename() const {
    return filename;
  }

  CountOutputFormat get_output_format() const {
    return output_format;
  }

  const std::string& get_column_name(const uint column_index) const {
    assert(column_index < column_names.size());
    return column_names[column_index];
  }

  // open file, return false if file could not be opened and error_is_fatal is false
  bool open(bool error_is_fatal = true);

  // flush buffer, open output file if needed, keep file open afterwards
  void flush();

private:
  void write_gdat_header();

  CountOutputFormat output_format;

  // name of the output file with the full path
  std::string filename;

  // does not include 'time' column
  // unused when output_format is CountOutputFormat::DAT
  std::vector<std::string> column_names;

  // number of rows to be stored, automatically flushes afterwards
  size_t buffer_size;

  // output stream
  std::ofstream fout;

  // buffer columns, size of this vector is the
  // same as column_names size
  std::vector<CountItemVector> data;

  bool open_for_append;
};

typedef std::vector<CountBuffer> CountBufferVector;

} /* namespace MCell */

#endif /* SRC4_COUNT_BUFFER_H_ */
