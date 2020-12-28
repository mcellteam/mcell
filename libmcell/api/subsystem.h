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

#ifndef API_SUBSYSTEM_H
#define API_SUBSYSTEM_H

#include <functional>

#include "generated/gen_subsystem.h"
#include "api/common.h"
#include "api/api_utils.h"
#include "api/species.h"
#include "api/surface_class.h"
#include "api/reaction_rule.h"

namespace BNG {
class ElemMolType;
class RxnRule;
class Cplx;
}

namespace MCell {
namespace API {

class Complex;

class Subsystem: public GenSubsystem {
public:
  SUBSYSTEM_CTOR()

  // from generated template
  void add_species(std::shared_ptr<Species> s) override {
    append_to_vec_canonical_name(species, s);
  }

  std::shared_ptr<Species> find_species(const std::string& name) override {
    return vec_find_by_name(species, name);
  }

  void add_reaction_rule(std::shared_ptr<ReactionRule> r) override {
    // reactions don't have to have name,
    // not checking whether same unnamed rule already exists
    append_to_vec_canonical_name(reaction_rules, r);
  }

  std::shared_ptr<ReactionRule> find_reaction_rule(const std::string& name) override {
    return vec_find_by_name(reaction_rules, name);
  }

  void add_surface_class(std::shared_ptr<SurfaceClass> sc) override {
    append_to_vec(surface_classes, sc);
  }

  std::shared_ptr<SurfaceClass> find_surface_class(const std::string& name) override {
    return vec_find_by_name(surface_classes, name);
  }

  void add_elementary_molecule_type(std::shared_ptr<ElementaryMoleculeType> mt) override {
    append_to_vec_canonical_name(elementary_molecule_types, mt);
  }

  std::shared_ptr<ElementaryMoleculeType> find_elementary_molecule_type(const std::string& name) override {
    return vec_find_by_name(elementary_molecule_types, name);
  }


  void load_bngl_molecule_types_and_reaction_rules(
      const std::string& file_name,
      const std::map<std::string, float_t>& parameter_overrides = std::map<std::string, float_t>()
  ) override;

  // added manually

  // go through all species and make sure that
  // 1) all identical elementary molecule type objects exist only once and
  // 2) all are also in the elementary_molecule_types array
  void unify_and_register_elementary_molecule_types();

  void dump() const;

  // auxiliary methods used also from Instantiation
  std::shared_ptr<API::Complex> convert_cplx(
      const BNG::BNGData& bng_data,
      const BNG::Cplx& bng_inst);

  std::shared_ptr<API::Complex> convert_cplx_w_orientation(
      const BNG::BNGData& bng_data,
      const BNG::Cplx& bng_inst,
      const Orientation orientation);

protected:
  void convert_bng_data_to_subsystem_data(const BNG::BNGData& bng_data);

private:
  void convert_elementary_molecule_type(const BNG::BNGData& bng_data, const BNG::ElemMolType& bng_mt);
  void convert_reaction_rule(const BNG::BNGData& bng_data, const BNG::RxnRule& bng_rr);
};

} // namespace API
} // namespace MCell

#endif // API_SUBSYSTEM_H
