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

#include "generated/gen_run_utils.h"

#include "bng/filesystem_include.h"

#include "api_utils.h"
#include "bng/bng_defines.h"
#include "src4/simulation_config.h"

using namespace std;

namespace MCell {
namespace API {

namespace run_utils {


std::string get_last_checkpoint_dir(const int seed) {
  // only supported parameter is seed for now

  string chkpt_seed_dir =
      DEFAULT_CHECKPOINTS_DIR + BNG::PATH_SEPARATOR + get_seed_dir_name(seed);

  string max_it_dir_name = "";

  // do we have a checkpoints/seed directory?
  if (fs::exists(chkpt_seed_dir)) {

    // find the highest iteration value
    int max = -1;
    for (const auto& entry :fs::directory_iterator(chkpt_seed_dir)) {
      string dir_name = entry.path().filename().string();

      // starts with "it_"?
      if (dir_name.find(DEFAULT_ITERATION_DIR_PREFIX) == 0) {
        string it_nr_str = dir_name.substr(
            DEFAULT_ITERATION_DIR_PREFIX.size(), dir_name.size() - DEFAULT_ITERATION_DIR_PREFIX.size());

        int it_nr;
        try {
          it_nr = stoi(it_nr_str);
        }
        catch (invalid_argument&) {
          continue;
        }

        if (it_nr > max) {
          // remember the directory with the highest iteration number so far
          max = it_nr;
          max_it_dir_name = entry.path().filename().string();
        }
      }
    }
  }

  if (max_it_dir_name != "") {
    return chkpt_seed_dir + BNG::PATH_SEPARATOR + max_it_dir_name + BNG::PATH_SEPARATOR;
  }
  else {
    return "";
  }
}


std::vector<std::string> remove_cwd(const std::vector<std::string> paths) {

  std::vector<std::string> res;

  // skip everything that can be interpreted as current directory
  string cwd = fs::current_path().string();
  for (const string& p: paths) {
    if (p != "" && p != "." && p != cwd) {
      res.push_back(p);
    }
  }

  return res;
}

} // namespace run_utils

} // namespace API
} // namespace MCell
