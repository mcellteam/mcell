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

#include "count_buffer.h"

#include "logging.h"

namespace MCell {

void CountItem::write(std::ostream& out) const {
  out << time << " " << value << "\n";
}

void CountBuffer::flush() {
  if (!fout.is_open()) {
    open(true);
  }

  for (const auto& item: data) {
    item.write(fout);
  }
  data.clear();
}


bool CountBuffer::open(bool error_is_fatal) {

  if (::make_parent_dir(filename.c_str()) ) {
    mcell_error("Failed to create parent directory for molecule count output.");
  }

  if (!open_for_append) {
    // create an empty file so that we know that nothing was stored
    fout.open(filename);
  }
  else {
    // appending is used when restoring a checkpoint
    // opens a new file if the file does not exist
    fout.open(filename, std::ofstream::out | std::ofstream::app);
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
