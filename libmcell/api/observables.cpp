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


CountOutputFormat Observables::count_output_format_from_path_or_file(const std::string& path_or_file) {
  const string gdat = ".gdat";
  size_t gdat_sz = gdat.size();
  size_t sz = path_or_file.size();
  if (sz > gdat.size() && path_or_file.substr(sz - gdat_sz) == gdat) {
    return CountOutputFormat::GDAT;
  }
  else {
    return CountOutputFormat::DAT;
  }
}


CountOutputFormat Observables::get_count_output_format(
    const CountOutputFormat observables_output_format,
    const std::string& observables_path_or_file
) {

  CountOutputFormat fmt;
  if (observables_output_format == CountOutputFormat::AUTOMATIC_FROM_EXTENSION) {
    return count_output_format_from_path_or_file(observables_path_or_file);
  }
  else {
    if (observables_path_or_file == "" &&
        observables_output_format == CountOutputFormat::GDAT) {
      throw RuntimeError(S("Attribute ") + NAME_OBSERVABLES_PATH_OR_FILE + " must not be empty " +
          "when " + NAME_OBSERVABLES_OUTPUT_FORMAT + " is " + NAME_ENUM_COUNT_OUTPUT_FORMAT + "." +
          NAME_EV_GDAT + ".");
    }
    return observables_output_format;
  }
}


void Observables::load_bngl_observables(
    const std::string& file_name,
    const std::string& observables_path_or_file,
    const std::map<std::string, double>& parameter_overrides,
    const CountOutputFormat observables_output_format) {
  BNG::BNGData bng_data;

  int num_errors = BNG::parse_bngl_file(file_name, bng_data, parameter_overrides);
  if (num_errors != 0) {
    throw RuntimeError("Could not parse BNGL file " + file_name + ".");
  }


  // now convert everything we parsed into the API classes so that the user can
  // inspect or manipulate it if needed
  convert_bng_data_to_observables_data(
      bng_data,
      get_count_output_format(observables_output_format, observables_path_or_file),
      observables_path_or_file);
}


void Observables::convert_bng_data_to_observables_data(
    const BNG::BNGData& bng_data,
    const CountOutputFormat observables_output_format,
    const std::string& observables_path_or_file) {

  for (const Observable& o: bng_data.get_observables()) {
    convert_observable(o, bng_data, observables_path_or_file, observables_output_format);
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
    const std::string& observables_path_or_file,
    const CountOutputFormat observables_output_format) {

  shared_ptr<API::Count> count = make_shared<Count>(DefaultCtorArgType());
  count->expression = make_shared<CountTerm>(DefaultCtorArgType());

  count->name = o.name;
  if (observables_output_format == CountOutputFormat::DAT) {
    if (is_set(observables_path_or_file)) {
      count->file_name = observables_path_or_file + o.name + ".dat";
    }
    else {
      // will be set during conversion
      count->file_name = STR_UNSET;
    }

  }
  else if (observables_output_format == CountOutputFormat::GDAT) {
    count->file_name = observables_path_or_file;
  }
  else {
    release_assert(false && "Invalid count output format.");
  }

  count->output_format = observables_output_format;

  release_assert(!o.patterns.empty() && "Observable has no pattern");
  if (o.patterns.size() == 1) {
    set_count_molecules_or_species_pattern(o.type, o.patterns[0], bng_data, *count->expression);
  }
  else {
    shared_ptr<API::CountTerm> top_level_term = nullptr;
    for (const auto& pat: o.patterns) {
      shared_ptr<API::CountTerm> term = make_shared<CountTerm>(DefaultCtorArgType());

      set_count_molecules_or_species_pattern(o.type, pat, bng_data, *term);

      if (is_set(top_level_term)) {
        top_level_term = top_level_term->__add__(term);
      }
      else {
        top_level_term = term;
      }
    }
    count->expression = top_level_term;
  }

  // remember for conversion in mcell4 coverter
  counts.push_back(count);
}


} // namespace API
} // namespace MCell

