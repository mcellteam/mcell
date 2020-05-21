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

#ifndef API_OBSERVABLES_H
#define API_OBSERVABLES_H

#include "generated/gen_observables.h"
#include "api/common.h"
#include "api/api_utils.h"
#include "api/viz_output.h"
#include "api/count.h"

namespace MCell {
namespace API {

class Observables: public GenObservables {
public:

  void add_viz_output(std::shared_ptr<VizOutput> viz_output) override {
    append_to_vec(viz_outputs, viz_output);
  };

  void add_count(std::shared_ptr<Count> count) override {
    append_to_vec(counts, count);
  };

  // added manually
  void dump() const {}
};

} // namespace API
} // namespace MCell

#endif // API_OBSERVABLES_H