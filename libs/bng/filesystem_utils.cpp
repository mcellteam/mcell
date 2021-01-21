/*
 * filesystem_utils.cpp
 *
 *  Created on: Jan 21, 2021
 *      Author: Adam
 */

#include "filesystem_utils.h"

#include <iostream>
#include <cassert>

#include "filesystem_include.h"

using namespace std;

void make_dir_for_file_w_multiple_attempts(const std::string& file_path) {
  string dir_path = fs::path(file_path).parent_path().string();

  make_dir_w_multiple_attempts(dir_path);
}


void make_dir_w_multiple_attempts(const std::string& dir_path) {
  if (fs::is_directory(dir_path)) {
    return;
  }

  // create parent dir if the parent does not exist
  fs::path parent_path = fs::path(dir_path).parent_path();
  if (!parent_path.empty() && !fs::is_directory(parent_path)) {
    make_dir_w_multiple_attempts(parent_path.string());
  }

  // parallel runs may collide when creating the directories,
  // trying it multiple times
  int num_attemts = 0;
  bool ok;
  do {
    ok = fs::create_directory(dir_path);

    if (!ok) {
      cerr << "Could not create directory '" << dir_path << "', trying again after 1s.\n";
      num_attemts++;
      errno = 0;
#ifndef _MSC_VER
      sleep(1);
#else
      _sleep(1);
#endif
      // another process might have created the directory
      ok = fs::is_directory(dir_path);
    }
  } while (!ok && num_attemts < 3);

  if (!ok) {
    cerr << "Could not create directory '" << dir_path << "', terminating simulation.\n";
    exit(1);
  }
}
