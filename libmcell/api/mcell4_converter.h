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

#ifndef LIBMCELL_API_MCELL4_CONVERTER_H_
#define LIBMCELL_API_MCELL4_CONVERTER_H_

#include <vector>
#include "defines.h"
#include "mol_or_rxn_count_event.h"

namespace BNG {
class ComponentType;
class Component;
class ElemMolType;
class ElemMol;
class Cplx;
class Species;
}

namespace MCell {

class World;
class Partition;
class InitialRegionMolecules;
class GeometryObject;
class RegionExprNode;
class ReleaseEvent;
class MolOrRxnCountTerm;
class MolOrRxnCountItem;
class Region;

namespace API {


class Model;
class ComponentType;
class Component;
class ElementaryMoleculeType;
class ElementaryMolecule;
class Complex;
class SurfaceProperty;
class InitialSurfaceRelease;
class Region;
class SurfaceRegion;
class GeometryObject;
class MoleculeReleaseInfo;
class ReleaseSite;
class CountTerm;
class Species;
class SurfaceClass;
class RngState;

class MCell4Converter {
public:
  MCell4Converter(Model* model_, World* world_) :
    model(model_), world(world_){
  }

  // convert all items in Model into the World representation
  // throws exception if anything went wrong
  // modifies model as well where it stores information for cases
  // when a value such a reaction rate was updated by the user
  void convert();

  // converter can be also used to convert individual objects
  // throws exception if anything went wrong
  MCell::ReleaseEvent* convert_single_release_event(
      const std::shared_ptr<API::ReleaseSite>& r);

  static void convert_rng_state(std::shared_ptr<RngState>& src, rng_state& dst);

private:
  species_id_t get_species_id(API::Species& s, const std::string class_name, const std::string object_name);
  species_id_t get_species_id_for_complex(API::Complex& ci, const std::string error_msg, const bool check_orientation = true);

  void get_geometry_bounding_box(Vec3& llf, Vec3& urb);
  void convert_simulation_setup();

  void convert_elementary_molecule_types();
  void convert_species();

  void convert_surface_class_rxn(API::SurfaceProperty& sp, const BNG::Species& surface_reactant);
  void convert_surface_classes();

  BNG::component_type_id_t convert_component_type(API::ComponentType& api_ct);
  BNG::elem_mol_type_id_t convert_elementary_molecule_type(API::ElementaryMoleculeType& mt, const bool in_rxn_or_observables = false);

  BNG::Component convert_component_instance(API::Component& api_ci);
  BNG::ElemMol convert_molecule_instance(API::ElementaryMolecule& mi, const bool in_rxn_or_observables = false);
  BNG::Cplx convert_complex(API::Complex& inst, const bool in_observables = false, const bool in_rxn = false);

  void check_surface_compartments(const BNG::RxnRule& r, BNG::compartment_id_t& surf_comp_id);
  void set_vol_rxn_substance_orientation_from_compartment(
      BNG::RxnRule& r, const BNG::Compartment& surf_comp, BNG::Cplx& substance);
  void check_intermembrane_surface_reaction(const BNG::RxnRule& rxn);
  void convert_rxns();

  MCell::wall_index_t convert_wall_and_add_to_geom_object(
      const API::GeometryObject& src_obj, const uint side,
      MCell::Partition& p, MCell::GeometryObject& dst_obj
  );

  void convert_initial_surface_releases(
      const std::vector<std::shared_ptr<API::InitialSurfaceRelease>>& api_releases,
      std::vector<MCell::InitialRegionMolecules>& mcell_releases
  );

  void convert_concentration_clamp_release(
      const partition_id_t partition_id, const API::SurfaceClass& surface_class, const MCell::Region& mcell_region);

  MCell::region_index_t convert_surface_region(
      MCell::Partition& p,
      API::SurfaceRegion& surface_region, API::GeometryObject& o,
      MCell::GeometryObject& obj
  );
  void convert_geometry_objects();

  void check_surface_compartment_name_collision(const std::string& surface_compartment_name);
  void convert_compartments();

  MCell::RegionExprNode* convert_region_expr_recursively(
      const std::shared_ptr<API::Region>& region,
      MCell::ReleaseEvent* rel_event
  );
  void convert_region_expr(API::ReleaseSite& rel_site, MCell::ReleaseEvent* rel_event);
  void convert_molecule_list(
      const std::vector<std::shared_ptr<MoleculeReleaseInfo>>& molecule_list,
      const std::string& rel_site_name,
      MCell::ReleaseEvent* rel_event);
  void convert_release_events();

  void convert_count_term_leaf_and_init_counting_flags(
      const std::shared_ptr<API::CountTerm> ct, const int sign,
      MolOrRxnCountTermVector& terms
  );
  void convert_count_terms_recursively(
      const std::shared_ptr<API::CountTerm> ct,
      const int sign,
      MCell::MolOrRxnCountItem& info
  );
  void convert_mol_or_rxn_count_events_and_init_counting_flags();
  void init_species_counting_flags();

  void convert_viz_output_events();

  void convert_simulation_state();

  void add_ctrl_c_termination_event();

  void check_all_mol_types_have_diffusion_const();

  Model* model;
  World* world;
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_MCELL4_CONVERTER_H_ */
