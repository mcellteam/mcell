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

#ifndef API_OBSERVABLES_H
#define API_OBSERVABLES_H

#include "generated/gen_observables.h"
#include "api/api_common.h"
#include "api/api_utils.h"
#include "api/viz_output.h"
#include "api/count.h"

namespace BNG {
class Observable;
class BNGData;
}

namespace MCell {
namespace API {

class Observables: public GenObservables {
public:
  OBSERVABLES_CTOR()

  // TODO: add checks that each observable writes to a different file

  void add_viz_output(std::shared_ptr<VizOutput> viz_output) override {
    append_to_vec(viz_outputs, viz_output, true);
  };

  void add_count(std::shared_ptr<Count> count) override {
    append_to_vec(counts, count, true);
  };

  // throws error if viz_output is not present or if its output_files_prefix is not set
  // uses argument method_name for this error report
  std::string get_first_viz_output_files_prefix(const char* method_name);

  std::shared_ptr<Count> find_count(const std::string& name) override {
    return vec_find_by_name(counts, name);
  }

  void load_bngl_observables(
      const std::string& file_name,
      const std::string& observables_path_or_file = "",
      const std::map<std::string, double>& parameter_overrides = std::map<std::string, double>(),
      const CountOutputFormat observables_output_format = CountOutputFormat::AUTOMATIC_FROM_EXTENSION
  ) override;

  // added manually
  void dump() const;

protected:
  void convert_bng_data_to_observables_data(
      const BNG::BNGData& bng_data,
      const CountOutputFormat observables_output_format,
      const std::string& output_files_prefix
  );

  CountOutputFormat get_count_output_format(
      const CountOutputFormat observables_output_format,
      const std::string& observables_path_or_file
  );
private:
  void convert_observable(
      const BNG::Observable& o,
      const BNG::BNGData& bng_data,
      const std::string& observables_path_or_file,
      const CountOutputFormat observables_output_format
  );

  CountOutputFormat count_output_format_from_path_or_file(const std::string& path_or_file);
};

} // namespace API
} // namespace MCell

#endif // API_OBSERVABLES_H
