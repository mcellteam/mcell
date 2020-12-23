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

#include "api/python_exporter.h"
#include "api/python_export_constants.h"
#include "api/python_export_utils.h"
#include "api/model.h"
#include "world.h"

namespace MCell {
namespace API {

PythonExporter::PythonExporter(Model* model_) :
  model(model_) {

  assert(model != nullptr);
  world = model->get_world();
  assert(world != nullptr);
}


void PythonExporter::save_checkpoint(const std::string& dir) {
  // parameters
  // - includes rng state

  // subsystem


  // geometry

  // instantiation

  // observables

  // molecules
  // - volume
  // - surface

  // model
  // - config
  // - warnings
  // - notifications
}


} // namespace API
} // namespace MCell

