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

#ifndef API_COUNT_H
#define API_COUNT_H

#include "generated/gen_count.h"
#include "api/api_common.h"
#include "api/count_term.h"


namespace MCell {

class MolOrRxnCountEvent;

namespace API {

class Count: public GenCount {
public:
  COUNT_CTOR()

  void set_all_custom_attributes_to_default() override {
    count_event = nullptr;
  }

  void postprocess_in_ctor() override;
  void check_semantics() const override;

  double get_current_value() override;

  // count event, owned by Scheduler if every_n_timesteps > 0,
  // owned by World if every_n_timesteps == 0
  MolOrRxnCountEvent* count_event;

private:
  void set_automatically_output_format_if_needed();
};

} // namespace API
} // namespace MCell

#endif // API_COUNT_H
