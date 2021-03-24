
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h> // mkdir

#include "bng/shared_defines.h"
#include "bng/bng_defines.h" // includes bng_config.h
#include "bng/filesystem_utils.h"

using namespace std;

namespace BNG {

void BNGNotifications::dump() const {
  cout << "BNGNotifications:\n";
  cout << "  rxn_probability_changed: \t\t" << rxn_probability_changed << "\n";
}


float_t BNGConfig::get_default_rx_radius_3d() const {
  return (1.0 / sqrt_f(BNG_PI * grid_density)) / length_unit;
}


void BNGConfig::dump() const {
  cout << "BNGConfig:\n";
  cout << "  initial_seed: \t\t" << initial_seed << "\n";
  cout << "  time_unit: \t\t" << time_unit << "\n";
  cout << "  length_unit: \t\t" << length_unit << "\n";
  cout << "  grid_density: \t\t" << grid_density << "\n";
  cout << "  rx_radius_3d: \t\t" << rx_radius_3d << "\n";
  cout << "  rxn_and_species_report: \t\t" << rxn_and_species_report << "\n";
  // TODO: add dumps for BNGNotificatiosn and BNGWarnings

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


void BNGConfig::initialize_bng_report_files() {
  remove_report_file(get_warnings_report_file_name(), "Warnings");
  if (rxn_and_species_report) {
    initialize_report_file(get_rxn_report_file_name(), "RXN");
    initialize_report_file(get_species_report_file_name(), "Species");
  }
}


void BNGConfig::print_final_warnings() const {
  if (warnings.bimol_rxn_probability_over_1) {
    cerr <<
        "Warning: There was a bimolecular reaction with probability p > 1.0. This means that some reactions were missed. " <<
        "See report file " << get_warnings_report_file_name() << " for more details. A shorter time step may be needed.\n" <<
        "Additional details may be found in " << get_rxn_report_file_name() <<
        ", if it does not exist, it can be enabled by setting 'rxn_and_species_report' to true.\n";

  }
  if (warnings.warn_on_bimol_rxn_probability_over_05_less_1 &&
      warnings.bimol_rxn_probability_over_05_less_1) {
    cerr <<
        "Warning: There was a bimolecular reaction with probability p > 0.5 and p < 1.0. " <<
        "For best results, the probability of bimolecular reactions should be <= 0.5. " <<
        "See report file " << get_rxn_report_file_name() <<
        ", if it does not exist, it can be enabled by setting 'rxn_and_species_report' to true.\n";
  }
}


// get current date/time, format is YYYY-MM-DD HH:mm:ss
std::string get_current_date_time() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}


void append_to_report(const std::string& report_fname, const std::string& msg) {
  // check whether the directory exists and make it if not
  FSUtils::make_dir_for_file_w_multiple_attempts(report_fname);

  // then open or append
  ofstream of;
  of.open(report_fname, fstream::out | fstream::app);
  // not printing warning when file count not be opened
  if (of.is_open()) {
    of << msg;
    of.close();
  }
}

void initialize_report_file(const std::string& fname, const char* report_name) {
  // first check whether the directory exists and make it
  FSUtils::make_dir_for_file_w_multiple_attempts(fname);

  ofstream of;
  of.open(fname, fstream::out);
  if (of.is_open()) {
    of << report_name << " report, " << get_current_date_time() << "\n\n";
    of.close();
  }
  else {
    cout << "Could not open file " << fname << " for report generation, ignored.\n";
  }
}


void remove_report_file(const std::string& fname, const char* report_name) {
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


} // namespace BNG
