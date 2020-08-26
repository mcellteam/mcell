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

#include "nfsim_species_unifier.h"

#include <iostream>
#include <fstream>

using namespace std;

extern int bngldebug;

namespace MCell {

bool NFSimSpeciesUnifier::read_species_file(const std::string& input_file) {

  //bngldebug = 1;

  ifstream ifs;
  ifs.open(input_file, ifstream::in);
  if (!ifs.is_open()) {
    cerr << "Could not open input file " << input_file << ".\n";
    return false;
  }

  string line;
  uint linenr = 0;
  while (getline(ifs, line)) {
    linenr++;
    if (line.empty()) {
      continue;
    }
    if (line[0] == '#') {
      // comment
      continue;
    }

    // expecting that each line will be in this form:
    // BNGL_COMPLEX count
    size_t cplx_end = line.find(' ');
    if (cplx_end == string::npos) {
      cerr << "Did not find a separator ' ' on line " << linenr << ".\n";
      return false;
    }

    string cplx_string = line.substr(0, cplx_end);
    string count_string = line.substr(cplx_end);

    // parse cplx string
    BNG::CplxInstance cplx_inst(&bng_engine.get_data());
    int num_errors = BNG::parse_single_cplx_instance_string(
        cplx_string, bng_engine.get_data(),
        cplx_inst
    );
    if (num_errors != 0) {
      cerr << "Could not parse complex instance '" << cplx_string << "' on line " << linenr << ".\n";
      ifs.close();
      return false;
    }
    assert(!cplx_inst.mol_instances.empty());

    // see what species it is
    species_id_t id = bng_engine.get_all_species().find_or_add(
        BNG::Species(cplx_inst, bng_engine.get_data(), bng_engine.get_config(), false)
    );

    // convert count
    int count;
    try {
      count = stoi(count_string);
    }
    catch (const std::invalid_argument& ia) {
      cerr << "Could not convert count '" << count_string << "' on line " << linenr << ".\n";
      ifs.close();
      return false;
    }

    // remember this count
    auto it = counts_per_unique_species.find(id);
    if (it == counts_per_unique_species.end()) {
      counts_per_unique_species[id] = count;
    }
    else {
      counts_per_unique_species[id] += count;
    }
  }

  return true;
}


void NFSimSpeciesUnifier::print_unified_species() {
  for (auto it: counts_per_unique_species) {
    BNG::Species& s = bng_engine.get_all_species().get(it.first);
    cout << s.to_str(bng_engine.get_data()) << " " << it.second << "\n";
  }
}

} // namespace MCell