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

#ifndef API_GEN_SPECIES_H
#define API_GEN_SPECIES_H

#include "../api/common.h"
#include "../api/complex_instance.h"


namespace MCell {
namespace API {

class MoleculeInstance;

#define SPECIES_CTOR() \
    Species( \
        const std::string& name_, \
        const std::vector<std::shared_ptr<MoleculeInstance>> molecule_types_ = std::vector<std::shared_ptr<MoleculeInstance>>() \
    )  : GenSpecies(molecule_types_) { \
      class_name = "Species"; \
      name = name_; \
      molecule_types = molecule_types_; \
    }

class GenSpecies: public ComplexInstance {
public:
  GenSpecies( 
      const std::vector<std::shared_ptr<MoleculeInstance>> molecule_types_ = std::vector<std::shared_ptr<MoleculeInstance>>() 
  )  : ComplexInstance(molecule_types_)  {
  }
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  // --- methods ---
}; // GenSpecies

class Species;
py::class_<Species> define_pybinding_Species(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SPECIES_H
