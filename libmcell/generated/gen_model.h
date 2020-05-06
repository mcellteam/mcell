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

#ifndef API_GEN_MODEL_H
#define API_GEN_MODEL_H

#include "../api/common.h"

namespace MCell {
namespace API {

class GeometryObject;
class InstantiationData;
class ReactionRule;
class ReactionRule;
class ReleaseSite;
class ReleaseSite;
class Species;
class Species;
class Subsystem;

class GenModel {
public:
  virtual ~GenModel() {}
  // --- attributes ---
  // --- methods ---
  virtual void run_iterations(const long iterations) = 0;
  virtual void add_subsystem(std::shared_ptr<Subsystem> subsystem) = 0;
  virtual void add_instantiation_data(std::shared_ptr<InstantiationData> instantiation_data) = 0;
}; // GenModel

class Model;
py::class_<Model> define_pybinding_Model(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MODEL_H
