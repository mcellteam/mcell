/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef API_MODEL_H
#define API_MODEL_H

#include "../generated/gen_model.h"
#include "../api/common.h"
#include "subsystem.h"
#include "instantiation_data.h"
#include "config.h"
#include "warnings.h"
#include "notifications.h"

namespace MCell {
namespace API {

class Model: public GenModel, public Subsystem, public InstantiationData {
public:

  // from generated template
  void run_iterations(const long iterations) override {}
  void add_subsystem(std::shared_ptr<Subsystem> subsystem) override {}
  void add_instantiation_data(std::shared_ptr<InstantiationData> instantiation_data) override {}

  // added manually
  void dump() const;
};

} // namespace API
} // namespace MCell

#endif // API_MODEL_H
