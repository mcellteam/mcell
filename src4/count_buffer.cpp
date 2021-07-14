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

#include "count_buffer.h"

#include <iomanip>
#include <sstream>

#include "logging.h"
#include "util.h"

#include "bng/filesystem_utils.h"

using namespace std;

const uint GDAT_COLUMN_WIDTH = 14;


namespace MCell {

void CountItem::write_as_dat(std::ostream& out) const {
  out << time << " " << value << "\n";
}


static void write_gdat_value(std::ostream& outs, double d) {
  // there is no way to set number of digits in an exponent
  // so we must do it manually
  // split exponent and base as string - we do not want to do
  // any numerical computation here
  stringstream ss;
  ss << scientific << setprecision(GDAT_COLUMN_WIDTH - 6) << d;

  string s = ss.str();
  size_t pos = s.find('e');
  assert(pos != string::npos);

  string base = s.substr(0, pos);
  int exponent = stoi(s.substr(pos + 1));

  string sign;
  if (exponent >= 0) {
    sign = "+";
  }
  else {
    // setfill/setw does not add '0' for negative values
    sign = "-";
    exponent = -exponent;
  }

  outs << base << "e" << sign << setfill('0') << setw(2) << exponent;
}


void CountBuffer::flush() {
  if (!fout.is_open()) {
    open(true);
  }

  if (output_format == CountOutputFormat::DAT) {
    // there is a single column
    assert(data.size() == 1);
    for (const auto& item: data[0]) {
      assert(item.column_index == 0);
      item.write_as_dat(fout);
    }
  }
  else {
    assert(data.size() >= 1);

    // output row
    // expecting that each column has the same depth
    size_t num_rows = data[0].size();
    for (size_t row = 0; row < num_rows; row++) {

      // simply use time from the first column
      double time = data[0][row].time;
      fout << " ";
      write_gdat_value(fout, time);

      // output each column
      for (size_t col = 0; col < data.size(); col++) {
        assert(row < data[col].size());
        const auto& item = data[col][row];
        release_assert(cmp_eq(item.time, time, SQRT_EPS) && "Mismatch in gdat column times");
        fout << "  ";
        write_gdat_value(fout, item.value);
      }
      fout << "\n";
    }
  }

  fout.flush(); // flush the data so the user can see them
  for (auto& column: data) {
    column.clear();
  }
}


void CountBuffer::write_gdat_header() {
  assert(fout.is_open());
  assert(output_format == CountOutputFormat::GDAT);
  fout << "#";

  string time = "time";
  time.insert(time.begin(), GDAT_COLUMN_WIDTH - time.size(), ' ');
  fout << time;

  for (const string& name: column_names) {
    string col = name;
    if (GDAT_COLUMN_WIDTH > col.size()) {
      col.insert(col.begin(), GDAT_COLUMN_WIDTH - col.size() + 2, ' ');
    }
    else {
      col = " " + col;
    }
    fout << col;
  }
  fout << "\n";
}


bool CountBuffer::open(bool error_is_fatal) {

  FSUtils::make_dir_for_file_w_multiple_attempts(filename);

  if (!open_for_append) {
    // create an empty file so that we know that nothing was stored
    fout.open(filename);
  }
  else {
    // appending is used when restoring a checkpoint
    // opens a new file if the file does not exist
    fout.open(filename, std::ofstream::out | std::ofstream::app);
  }

  // write header
  if (output_format == CountOutputFormat::GDAT && !open_for_append) {
    write_gdat_header();
  }

  if (!fout.is_open()) {
    mcell_warn("Could not open file %s for writing.", filename.c_str());
    if (error_is_fatal) {
      mcell_error("Terminating due to error.");
    }
    return false;
  }
  return true;
}


void CountBuffer::flush_and_close() {

  flush();

  if (fout.is_open()) {
    fout.close();
  }
  else {
    bool ok = open(false);
    if (ok) {
      fout.close();
    }
  }
}

} /* namespace MCell */
