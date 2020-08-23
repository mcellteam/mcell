
#include <iostream>
#include <fstream>
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
  cout << "  debug_reactions: \t\t" << debug_reactions << " [bool] \t\t\n";

  notifications.dump();
}


std::string BNGConfig::get_rxn_report_file_name() const {
  return RXN_REPORT_PREFIX + to_string(initial_seed) + REPORT_EXT;
}


std::string BNGConfig::get_species_report_file_name() const {
  return SPECIES_REPORT_PREFIX + to_string(initial_seed) + REPORT_EXT;
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


void BNGConfig::initialize_report_files() {
  if (reporting) {

    ofstream of_rxn;
    of_rxn.open(get_rxn_report_file_name(), fstream::out);
    if (of_rxn.is_open()) {
      of_rxn << "RXN report, " << current_date_time() << "\n\n";
      of_rxn.close();
    }
    else {
      cout << "Could not open file " << get_rxn_report_file_name() << " for report generation, ignored.\n";
    }

    ofstream of_species;
    of_species.open(get_species_report_file_name(), fstream::out);
    if (of_species.is_open()) {
      of_species << "Species report, " << current_date_time() << "\n\n";
      of_species.close();
    }
    else {
      cout << "Could not open file " << get_species_report_file_name() << " for report generation, ignored.\n";
    }
  }
}


} // namespace BNG
