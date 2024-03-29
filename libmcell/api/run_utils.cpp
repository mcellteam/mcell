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

#include "generated/gen_run_utils.h"

#include "bng/filesystem_utils.h"

#include "api/api_utils.h"
#include "generated/gen_names.h"
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
  if (FSUtils::is_dir(chkpt_seed_dir)) {

    // find the highest iteration value
    int max = -1;
    std::vector<std::string> dirs;
    FSUtils::list_dir(chkpt_seed_dir, dirs);
    for (const string& dir_name: dirs) {

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
          max_it_dir_name = dir_name;
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
  string cwd = FSUtils::get_current_dir();
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
