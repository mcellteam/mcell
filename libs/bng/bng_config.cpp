
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>

#include "defines_shared.h" // includes also bng_config.h
#include "bng_defines.h" // includes also bng_config.h

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
  return RXN_REPORT_PREFIX + seed_as_str() + REPORT_EXT;
}


std::string BNGConfig::get_species_report_file_name() const {
  return SPECIES_REPORT_PREFIX + seed_as_str() + REPORT_EXT;
}


std::string BNGConfig::get_warnings_report_file_name() const {
  return WARNINGS_REPORT_PREFIX + seed_as_str() + REPORT_EXT;
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


static void initialize_file(const std::string& fname, const char* report_name) {
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


void BNGConfig::initialize_report_files() {
  initialize_file(get_warnings_report_file_name(), "Warnings");
  if (rxn_and_species_report) {
    initialize_file(get_rxn_report_file_name(), "RXN");
    initialize_file(get_species_report_file_name(), "Species");
  }
}


void append_to_report(const std::string& report_fname, const std::string& msg) {
  ofstream of;
  of.open(report_fname, fstream::out | fstream::app);
  // not printing warning when file count not be opened
  if (of.is_open()) {
    of << msg;
    of.close();
  }
}

} // namespace BNG
