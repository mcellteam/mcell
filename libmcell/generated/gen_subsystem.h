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

#ifndef API_GEN_SUBSYSTEM_H
#define API_GEN_SUBSYSTEM_H

#include "../api/mcell.h"

namespace MCell {
namespace API {

class Species;

class GenSubsystem {
public:
  virtual ~GenSubsystem() {}
  // --- attributes ---
  // --- methods ---
  virtual void add_species(const Species* s) = 0;
  virtual Species* find_species(const std::string& name) = 0;
}; // GenSubsystem

class Subsystem;
py::class_<Subsystem> define_pybinding_Subsystem(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SUBSYSTEM_H
