/*
 * filesystem_utils.cpp
 *
 *  Created on: Jan 21, 2021
 *      Author: Adam
 */

#include "filesystem_utils.h"
#include "bng_defines.h"

#include <iostream>
#include <cassert>

// when enabled, use posix implementations at all times
//#define POSIX_TEST

#if defined(__APPLE__) || defined(POSIX_TEST)

// apple supports std::filesystem from version 10.15
#define HAS_FS 0
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

#elif __GNUC__ < 8 && !defined(_MSC_VER)

// gcc 6 & 7 have filesystem still under the experimental features
#define HAS_FS 1
#include <unistd.h>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#else

#define HAS_FS 1
#include <filesystem>
namespace fs = std::filesystem;

#endif


using namespace std;

namespace FSUtils {

bool is_dir(const std::string& dir_path) {
#if HAS_FS
  return fs::is_directory(dir_path);
#else
  struct stat sb;
  if (stat(dir_path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
    return true;
  }
  else {
    return false;
  }
#endif
}


static bool make_dir(const std::string& dir_path) {
#if HAS_FS
  return fs::create_directory(dir_path);
#else
  return mkdir(dir_path.c_str(), 0777) == 0;
#endif
}


static std::string get_parent_path(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  if (pos != string::npos) {
    return path.substr(0, pos);
  }
  else {
    return "";
  }
}


void make_dir_for_file_w_multiple_attempts(const std::string& file_path) {
  string dir_path = get_parent_path(file_path);

  make_dir_w_multiple_attempts(dir_path);
}


static bool is_empty_path(const std::string& dir_path) {
  return dir_path == "" ||  dir_path == "." || dir_path == std::string(".") + BNG::PATH_SEPARATOR;
}


void make_dir_w_multiple_attempts(const std::string& dir_path) {

  if (is_empty_path(dir_path) || is_dir(dir_path)) {
    return;
  }

  // create parent dir if the parent does not exist
  string parent_path = get_parent_path(dir_path);
  if (!parent_path.empty() && !is_dir(parent_path)) {
    make_dir_w_multiple_attempts(parent_path);
  }

  // parallel runs may collide when creating the directories,
  // trying it multiple times
  int num_attemts = 0;
  bool ok;
  do {
    ok = make_dir(dir_path);

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
      ok = is_dir(dir_path);
    }
  } while (!ok && num_attemts < 3);

  if (!ok) {
    cerr << "Could not create directory '" << dir_path << "', terminating simulation.\n";
    exit(1);
  }
}


void list_dir(const std::string& dir_path, std::vector<std::string>& dirs) {
  dirs.clear();
#if HAS_FS
  for (const auto& entry :fs::directory_iterator(dir_path)) {
    string dir_name = entry.path().filename().string();
    dirs.push_back(dir_name);
  }
#else
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir (dir_path.c_str())) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir (dir)) != NULL) {
      dirs.push_back(ent->d_name);
    }
    closedir (dir);
  }
#endif
}


std::string get_current_dir() {
#if HAS_FS
   return fs::current_path().string();
#else
   char buf[2048];
   char* res = getcwd(buf, 2048);
   release_assert(res != nullptr && "Directory name is too long");
   string dir = buf;
   return dir;
#endif
}

} // namespace FSUtils

