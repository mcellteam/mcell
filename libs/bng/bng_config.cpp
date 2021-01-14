
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h> // mkdir

#include "bng/shared_defines.h"
#include "bng/bng_defines.h" // includes bng_config.h

using namespace std;

namespace BNG {

void BNGNotifications::dump() const {
  cout << "BNGNotifications:\n";
  cout << "  rxn_probability_changed: \t\t" << rxn_probability_changed << " [bool] \t\t\n";
}


void BNGConfig::dump() const {
  cout << "BNGConfig:\n";
  cout << "  initial_seed: \t\t" << initial_seed << " [uint] \t\t\n";
  cout << "  time_unit: \t\t" << time_unit << " [float_t] \t\t\n";
  cout << "  length_unit: \t\t" << length_unit << " [float_t] \t\t\n";
  cout << "  grid_density: \t\t" << grid_density << " [float_t] \t\t\n";
  cout << "  rx_radius_3d: \t\t" << rx_radius_3d << " [float_t] \t\t\n";
  cout << "  rxn_and_species_report: \t\t" << rxn_and_species_report << " [bool] \t\t\n";
  cout << "  bng_verbosity_level: \t\t" << bng_verbosity_level << " [bool] \t\t\n";

  notifications.dump();
}


std::string BNGConfig::seed_as_str() const {
  stringstream ss;
  ss << std::setw(5) << std::setfill('0') << initial_seed;
  return ss.str();
}

std::string BNGConfig::get_rxn_report_file_name() const {
  return std::string(REPORT_DIR) + PATH_SEPARATOR + RXN_REPORT_PREFIX + seed_as_str() + REPORT_EXT;
}


std::string BNGConfig::get_species_report_file_name() const {
  return std::string(REPORT_DIR) + PATH_SEPARATOR + SPECIES_REPORT_PREFIX + seed_as_str() + REPORT_EXT;
}


std::string BNGConfig::get_warnings_report_file_name() const {
  return std::string(REPORT_DIR) + PATH_SEPARATOR + WARNINGS_REPORT_PREFIX + seed_as_str() + REPORT_EXT;
}


// get current date/time, format is YYYY-MM-DD HH:mm:ss
static const std::string current_date_time() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}


static bool dir_exists(const string& path) {
  struct stat sb;
  if (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
    return true;
  }
  return false;
}


static void make_reports_dir(const string& file_path) {
  size_t pos = file_path.find_last_of("/\\");
  assert(pos != string::npos);
  string dir_path = file_path.substr(0, pos);

  if (dir_exists(dir_path)) {
    return;
  }

  // parallel runs may collide when creating the directories,
  // trying it multiple times
  int res;
  uint num_attemts = 0;
  do {
#ifndef _WIN64
    res = mkdir(dir_path.c_str(), 0777);
#else
    res = mkdir(dir_path.c_str());
#endif
    if (res != 0) {
      if (errno != EEXIST) {
        cerr << "Could not create directory '" << dir_path << "', trying again after 1s.\n";
        num_attemts++;
        errno = 0;
#ifndef _WIN64
        sleep(1);
#else
        // not sure, maybe gets optimized away
        int i = 0;
        while (i < 1000000) {
          i++;
        }
#endif
      }
      else {
        // directory now exists
        res = 0;
        errno = 0;
      }
    }
  } while (res != 0 && num_attemts < 3);

  if (res != 0) {
    cerr << "Could not create directory '" << dir_path << "', terminating simulation.\n";
    exit(1);
  }
}


static void initialize_file(const std::string& fname, const char* report_name) {
  // first check whether the directory exists and make it
  make_reports_dir(fname);

  ofstream of;
  of.open(fname, fstream::out);
  if (of.is_open()) {
    of << report_name << " report, " << current_date_time() << "\n\n";
    of.close();
  }
  else {
    cout << "Could not open file " << fname << " for report generation, ignored.\n";
  }
}


static void remove_file(const std::string& fname, const char* report_name) {
  ofstream of;
  of.open(fname, fstream::in);
  if (of.is_open()) {
    // file exists
    of.close();
    int res = remove(fname.c_str());
    if (res != 0){
      cout << "Could not remove file " << fname << " for report generation, ignored.\n";
    }
  }
}


void BNGConfig::initialize_report_files() {
  remove_file(get_warnings_report_file_name(), "Warnings");
  if (rxn_and_species_report) {
    initialize_file(get_rxn_report_file_name(), "RXN");
    initialize_file(get_species_report_file_name(), "Species");
  }
}


void append_to_report(const std::string& report_fname, const std::string& msg) {
  // check whether the directory exists and make it if not
  make_reports_dir(report_fname);

  // then open or append
  ofstream of;
  of.open(report_fname, fstream::out | fstream::app);
  // not printing warning when file count not be opened
  if (of.is_open()) {
    of << msg;
    of.close();
  }
}

} // namespace BNG
