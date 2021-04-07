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

#include "api/observables.h"
#include "api/subsystem.h"

#include "bng/bng.h"

#include "api/component_type.h"
#include "api/component.h"
#include "api/elementary_molecule_type.h"
#include "api/elementary_molecule.h"
#include "api/complex.h"

#include "api/api_utils.h"

using namespace std;
using namespace BNG;

namespace MCell {
namespace API {

void Observables::dump() const {
  std::cout << to_str() << "\n";
}


std::string Observables::get_first_viz_output_files_prefix(const char* method_name) {
  // use the first viz_output
  if (viz_outputs.empty()) {
    throw ValueError(
        S("Method ") + method_name + " requires a file argument when there are no instances of " +
        NAME_CLASS_VIZ_OUTPUT + " present."
    );
  }

  if (!is_set(viz_outputs[0]->output_files_prefix)) {
    throw ValueError(
        S("Method ") + method_name + ": the first VizOutput instance does not have its " +
        NAME_OUTPUT_FILES_PREFIX + " set."
    );
  }
  return viz_outputs[0]->output_files_prefix;
}


void Observables::load_bngl_observables(
    const std::string& file_name,
    const std::string& output_files_prefix,
    const std::map<std::string, float_t>& parameter_overrides) {
  BNG::BNGData bng_data;

  int num_errors = BNG::parse_bngl_file(file_name, bng_data, parameter_overrides);
  if (num_errors != 0) {
    throw RuntimeError("Could not parse BNGL file " + file_name + ".");
  }

  // now convert everything we parsed into the API classes so that the user can
  // inspect or manipulate it if needed
  convert_bng_data_to_observables_data(bng_data, output_files_prefix);
}


void Observables::convert_bng_data_to_observables_data(
    const BNG::BNGData& bng_data,
    const std::string& output_files_prefix) {

  for (const Observable& o: bng_data.get_observables()) {
    convert_observable(o, bng_data, output_files_prefix);
  }
}


static void set_count_molecules_or_species_pattern(
    const BNG::ObservableType type,
    const BNG::Cplx& bng_pattern,
    const BNG::BNGData& bng_data,
    API::CountTerm& count
) {

  std::shared_ptr<API::Complex> pattern =
      Complex::construct_from_bng_cplx(bng_data, bng_pattern);

  if (type == ObservableType::Molecules) {
    count.molecules_pattern = pattern;
  }
  else if (type == ObservableType::Species) {
    count.species_pattern = pattern;
  }
  else {
    release_assert(false);
  }
}


void Observables::convert_observable(
    const BNG::Observable& o,
    const BNG::BNGData& bng_data,
    const std::string& output_files_prefix) {

  shared_ptr<API::Count> count = make_shared<Count>(true);

  count->name = o.name;
  count->file_name = output_files_prefix + o.name + ".dat";

  release_assert(!o.patterns.empty() && "Observable has no pattern");
  if (o.patterns.size() == 1) {
    set_count_molecules_or_species_pattern(o.type, o.patterns[0], bng_data, *count);
  }
  else {
    shared_ptr<API::CountTerm> top_level_term = nullptr;
    for (const auto& pat: o.patterns) {
      shared_ptr<API::CountTerm> term = make_shared<CountTerm>();

      set_count_molecules_or_species_pattern(o.type, pat, bng_data, *term);

      if (is_set(top_level_term)) {
        top_level_term = top_level_term->__add__(term);
      }
      else {
        top_level_term = term;
      }
    }
    count->count_expression = top_level_term;
  }

  // remember for conversion in mcell4 coverter
  counts.push_back(count);
}


} // namespace API
} // namespace MCell

