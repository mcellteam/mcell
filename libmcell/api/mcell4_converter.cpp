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

#include "mcell4_converter.h"
#include "model.h"
#include "world.h"
#include "release_event.h"
#include "viz_output_event.h"
#include "mol_or_rxn_count_event.h"
#include "periodic_call_event.h"
#include "geometry.h"
#include "rng.h"
#include "isaac64.h"
#include "mcell_structs.h"
#include "bng/bng.h"

#include "api/mcell.h"
#include "api/bindings.h"

using namespace std;

namespace MCell {
namespace API {

static orientation_t convert_orientation(const Orientation o, const bool not_set_is_none = false) {
  switch (o) {
    case Orientation::DOWN:
      return ORIENTATION_DOWN;
    case Orientation::NONE:
      return ORIENTATION_NONE;
    case Orientation::UP:
      return ORIENTATION_UP;
    case Orientation::NOT_SET:
      if (not_set_is_none) {
        return ORIENTATION_NONE;
      }
      else {
        return ORIENTATION_NOT_SET;
      }
    case Orientation::ANY:
      return ORIENTATION_NONE;
    default:
      throw ValueError("Invalid Orientation value " + to_string((int)o) + ".");
  }
}


static viz_mode_t convert_viz_mode(const VizMode m) {
  switch (m) {
    case VizMode::ASCII:
      return ASCII_MODE;
    case VizMode::CELLBLENDER:
      return CELLBLENDER_MODE;
    default:
      throw ValueError("Invalid VizMode value " + to_string((int)m) + ".");
  }
}


void MCell4Converter::convert(Model* model_, World* world_) {
  model = model_;
  world = world_;

  convert_simulation_setup();
  convert_species();
  world->create_diffusion_events();

  convert_surface_classes();

  convert_rxns();
  init_species_rxn_flags();

  // at this point, we need to create the first (and for now the only) partition
  // create initial partition with center at 0,0,0
  partition_id_t index = world->add_partition(Vec3(0, 0, 0));
  assert(index == PARTITION_ID_INITIAL);

  convert_geometry_objects();

  // uses random generator state
  Geometry::check_for_overlapped_walls(world);

  convert_release_events();

  convert_mol_or_rxn_count_events_and_init_counting_flags();

  convert_viz_output_events();

  add_ctrl_c_termination_event();
}


float_t MCell4Converter::get_max_abs_dimension_of_any_vertex() {
  float_t max = 0;

  for (std::shared_ptr<API::GeometryObject>& o: model->geometry_objects) {
    // go through all vertices
    for (auto& v: o->vertex_list) {
      for (float_t dim: v) {
        float_t abs_dim = MCell::fabs_f(dim);
        if (abs_dim > max) {
          max = abs_dim;
        }
      }
    }
  }

  return max;
}


void MCell4Converter::convert_simulation_setup() {
  const API::Config& config = model->config;

  world->total_iterations = config.total_iterations_hint;
  world->config.time_unit = config.time_step;

  float_t grid_density = config.surface_grid_density;
  world->config.grid_density = grid_density;

  float_t length_unit = 1/sqrt_f(config.surface_grid_density);
  world->config.length_unit = length_unit;

  if (is_set(config.interaction_radius)) {
    // NOTE: mcell3 does not convert the unit of the interaction radius in parser which
    // seems a bit weird
    world->config.rx_radius_3d = config.interaction_radius / length_unit;
  }
  else {
    world->config.rx_radius_3d = (1.0 / sqrt_f(MY_PI * grid_density)) / length_unit;
  }

  float_t vacancy_search_dist = config.vacancy_search_distance / length_unit; // Convert units
  world->config.vacancy_search_dist2 = vacancy_search_dist * vacancy_search_dist; // and take square

  world->config.randomize_smol_pos = !config.center_molecules_on_grid;

  world->seed_seq = config.seed;
  rng_init(&world->rng, world->seed_seq);

  float_t max_dimension = get_max_abs_dimension_of_any_vertex();
  float_t auto_partition_dimension = (max_dimension + PARTITION_EDGE_EXTRA_MARGIN_UM) * 2;

  if (auto_partition_dimension > config.partition_dimension) {
    cout <<
        "Value of " << NAME_CLASS_MODEL << "." << NAME_CONFIG << "." << NAME_PARTITION_DIMENSION <<
        " " << config.partition_dimension <<
        " is smaller than the automatically determined value " << auto_partition_dimension <<
        ", using the automatic value.\n";
    world->config.partition_edge_length = auto_partition_dimension / length_unit;
  }
  else {
    world->config.partition_edge_length = config.partition_dimension / length_unit;
  }

  int num_subparts = config.partition_dimension / config.subpartition_dimension;
  assert(num_subparts > 0);
  if (num_subparts % 2 == 1) {
    // the number of subparts must be even
    num_subparts++;
  }
  world->config.num_subpartitions_per_partition = num_subparts;

  // this option in MCell3 was removed in MCell4
  world->config.use_expanded_list = true;

  // compute other constants
  world->config.init();
}


void MCell4Converter::convert_species() {
  for (std::shared_ptr<API::Species>& s: model->species) {
    BNG::Species new_species;
    new_species.name = s->name;

    bool is_vol = false;
    if (is_set(s->diffusion_constant_3d)) {
      is_vol = true;
      new_species.D = s->diffusion_constant_3d;
      new_species.set_is_vol();
      new_species.update_space_and_time_step(world->config.time_unit, world->config.length_unit);
    }
    else if (is_set(s->diffusion_constant_2d)) {
      is_vol = false;
      new_species.D = s->diffusion_constant_2d;
      new_species.set_is_surf();
      new_species.update_space_and_time_step(world->config.time_unit, world->config.length_unit);
    }
    else if (is_species_superclass(new_species.name)) {
      is_vol = new_species.name != ALL_SURFACE_MOLECULES;
      // these values are not really used, they are initialized for comparisons
      new_species.D = 0;
      new_species.space_step = 0;
      new_species.time_step = 0;
    }
    else
    {
      throw ValueError("Neither diffusion_constant_2d nor diffusion_constant_3d was set.");
    }
	
    new_species.set_flag(BNG::SPECIES_FLAG_CANT_INITIATE, s->target_only); // default is false

    // we must add a complex instance as the single molecule type in the new species
    // define a molecule type with no components
    BNG::MolType mol_type;
    mol_type.name = new_species.name; // name of the mol type is the same as for our species
    BNG::mol_type_id_t mol_type_id = world->bng_engine.get_data().find_or_add_molecule_type(mol_type);

    BNG::MolInstance mol_inst;
    mol_inst.mol_type_id = mol_type_id;
    if (is_vol) {
      mol_inst.set_is_vol();
    }
    else {
      mol_inst.set_is_surf();
    }

    new_species.mol_instances.push_back(mol_inst);

    new_species.finalize();
    species_id_t new_species_id = world->get_all_species().find_or_add(new_species);

    // remember which species we created
    s->species_id = new_species_id;

    // and also set superclass id for container if needed
    if (new_species.name == ALL_MOLECULES) {
      world->get_all_species().set_all_molecules_species_id(new_species_id);
    }
    else if (new_species.name == ALL_VOLUME_MOLECULES) {
      world->get_all_species().set_all_volume_molecules_species_id(new_species_id);
    }
    else if (new_species.name == ALL_SURFACE_MOLECULES) {
      world->get_all_species().set_all_surface_molecules_species_id(new_species_id);
    }

  }
}


void MCell4Converter::convert_surface_class_rxn(
    API::SurfaceProperty& sp, const BNG::Species& surface_reactant) {

  assert(sp.affected_species->species_id != SPECIES_ID_INVALID);
  BNG::Species& affected_species = world->get_all_species().get(sp.affected_species->species_id);

  BNG::RxnRule rxn;

  rxn.name = affected_species.name + "+" + surface_reactant.name;

  switch (sp.type) {
    case API::SurfacePropertyType::ABSORPTIVE:
      rxn.type = BNG::RxnType::Standard;
      break;
    case API::SurfacePropertyType::REFLECTIVE:
      rxn.type = BNG::RxnType::Reflect;
      break;
    case API::SurfacePropertyType::TRANSPARENT:
      rxn.type = BNG::RxnType::Transparent;
      break;
    default:
      throw ValueError("Invalid SurfaceProperty type for " + surface_reactant.name + ".");
  }

  // all these reactions happen always
  rxn.rate_constant = FLT_GIGANTIC;

  rxn.append_reactant(affected_species);
  rxn.append_reactant(surface_reactant); // copies the input reactant

  // this is the default orientation used for surfaces in a reaction
  rxn.reactants[1].set_orientation(ORIENTATION_UP);

  // add reaction and remember mapping
  sp.rxn_rule_id = world->get_all_rxns().add_finalized_no_update(rxn);
}


void MCell4Converter::convert_surface_classes() {
  for (std::shared_ptr<API::SurfaceClass>& sc: model->surface_classes) {
    // each surface class is represented by a species
    BNG::Species sc_species;
    sc_species.name = sc->name;

    sc_species.set_is_reactive_surface();
    // sets steps to 0
    sc_species.space_step = 0;
    sc_species.time_step = 0;
    sc_species.D = 0;

    // we must add a complex instance as the single molecule type in the new species
    // define a molecule type with no components
    BNG::MolType mol_type;
    mol_type.name = sc_species.name; // name of the mol type is the same as for our species
    BNG::mol_type_id_t mol_type_id = world->bng_engine.get_data().find_or_add_molecule_type(mol_type);

    BNG::MolInstance mol_inst;
    mol_inst.mol_type_id = mol_type_id;
    mol_inst.set_is_reactive_surface();
    sc_species.mol_instances.push_back(mol_inst);

    sc_species.finalize();
    species_id_t new_species_id = world->get_all_species().find_or_add(sc_species);

    // remember which species we created
    sc->species_id = new_species_id;

    // and we also need to add a reaction for each property
    if (sc->properties.empty()) {
      convert_surface_class_rxn(*dynamic_pointer_cast<SurfaceProperty>(sc), sc_species);
    }
    else {
      for (shared_ptr<SurfaceProperty>& sp: sc->properties) {
        convert_surface_class_rxn(*sp, sc_species);
      }
    }
  }
}


BNG::ComponentType MCell4Converter::convert_component_type(API::ComponentType& ct) {
  throw RuntimeError("Components are not supported yet");
}


BNG::ComponentInstance MCell4Converter::convert_component_instance(API::ComponentInstance& ci) {
  throw RuntimeError("Components are not supported yet");
}


BNG::MolType MCell4Converter::convert_molecule_type(API::ElementaryMoleculeType& mt) {
  BNG::MolType res;

  res.name = mt.name;
  if (!mt.components.empty()) {
    throw RuntimeError("Components are not supported yet");
  }

  return res;
}


BNG::MolInstance MCell4Converter::convert_molecule_instance(API::ElementaryMoleculeInstance& mi) {
  BNG::MolInstance res;

  BNG::MolType mt = convert_molecule_type(*mi.elementary_molecule_type);
  res.mol_type_id = world->bng_engine.get_data().find_or_add_molecule_type(mt);

  if (!mi.components.empty()) {
    throw RuntimeError("Components are not supported yet");
  }

  return res;
}


BNG::CplxInstance MCell4Converter::convert_complex_instance(API::ComplexInstance& inst) {
  // create a temporary cplx instance that we will use for search
  BNG::CplxInstance cplx_inst;

  for (std::shared_ptr<API::ElementaryMoleculeInstance>& m: inst.molecule_instances) {
    BNG::MolInstance mi = convert_molecule_instance(*m);

    cplx_inst.mol_instances.push_back(mi);
  }
  orientation_t orient = convert_orientation(inst.orientation, true);
  cplx_inst.set_orientation(orient);
  cplx_inst.finalize();

  // we need to find existing species that we match
  // at least for now until full BNG support will be ready
  species_id_t species_id = world->get_all_species().find_simple_species_id(cplx_inst);
  if (species_id == SPECIES_ID_INVALID) {
    throw ("Could not match reactant or product " + cplx_inst.to_str(world->bng_engine.get_data(), true) +
        " to any existing species.");
  }
  return world->bng_engine.create_cplx_instance_for_species(species_id, orient);
}


void MCell4Converter::convert_rxns() {
  for (std::shared_ptr<API::ReactionRule>& r: model->reaction_rules) {

    bool is_reversible = is_set(r->rev_rate);

    BNG::RxnRule rxn;

    if (is_set(r->name)) {
      rxn.name = r->name;
    }
    rxn.type = BNG::RxnType::Standard;

    if (is_set(r->fwd_rate)) {
      rxn.rate_constant = r->fwd_rate;
    }
    else {
      assert(!r->variable_rate.empty());
      for (auto& time_and_rate: r->variable_rate) {
        assert(time_and_rate.size() == 2);
        BNG::RxnRateInfo info;
        info.time = time_and_rate[0] / world->config.time_unit;
        info.rate_constant = time_and_rate[1];
        rxn.variable_rates.push_back(info);
      }
      rxn.update_variable_rxn_rate(0, nullptr);
    }

    for (std::shared_ptr<API::ComplexInstance>& rinst: r->reactants) {
      // convert to BNG::ComplexInstance using existing or new BNG::molecule_id

      BNG::CplxInstance reactant = convert_complex_instance(*rinst);
      rxn.append_reactant(reactant);
    }

    for (std::shared_ptr<API::ComplexInstance>& pinst: r->products) {
      // convert to BNG::ComplexInstance using existing or new BNG::molecule_id

      BNG::CplxInstance product = convert_complex_instance(*pinst);
      rxn.append_product(product);
    }

    // reverse reaction
    BNG::RxnRule rxn_rev;
    if (is_reversible) {
      rxn_rev.type = BNG::RxnType::Standard;
      if (is_set(r->rev_name)) {
        rxn.name = r->rev_name;
      }
      rxn_rev.rate_constant = r->rev_rate;
      rxn_rev.reactants = rxn.products;
      rxn_rev.products = rxn.reactants;
    }

    // add reaction(s) and remember mapping
    r->fwd_rxn_rule_id = world->get_all_rxns().add_finalized_no_update(rxn);

    if (is_reversible) {
      r->rev_rxn_rule_id = world->get_all_rxns().add_finalized_no_update(rxn_rev);
    }
  }
}



void MCell4Converter::init_species_rxn_flags() {
  BNG::SpeciesContainer& all_species = world->get_all_species();
  BNG::RxnContainer& all_rxns = world->get_all_rxns();


  species_id_t all_molecules_species_id = all_species.get_all_molecules_species_id();
  species_id_t all_volume_molecules_species_id = all_species.get_all_volume_molecules_species_id();
  species_id_t all_surface_molecules_species_id = all_species.get_all_surface_molecules_species_id();

  bool has_vol_vol_rxn = false;

  bool all_vol_mols_can_react_with_surface = false;
  bool all_surf_mols_can_react_with_surface = false;

  auto& species_vector = all_species.get_species_vector();

  // the first three classes are the superclasses
  assert(species_vector.size() > 3 &&
      species_vector[0].id == all_molecules_species_id &&
      species_vector[1].id == all_volume_molecules_species_id &&
      species_vector[2].id == all_surface_molecules_species_id
  );

  // set species flags (we already have all the reactions)
  for (BNG::Species& sp: species_vector) {
    // setup for ordinary molecules run after ALL_MOLECULES and ALL_VOLUME_MOLECULES
    // were processed
    if (sp.is_vol() && all_vol_mols_can_react_with_surface) {
      sp.set_flag(BNG::SPECIES_FLAG_CAN_VOLWALL);
    }

    if (sp.is_surf() && all_surf_mols_can_react_with_surface) {
      sp.set_flag(BNG::SPECIES_FLAG_CAN_REGION_BORDER);
    }

    // get reactions, this also creates all reaction classes for the species that we currently have
    BNG::SpeciesRxnClassesMap* rxn_classes =
        all_rxns.get_bimol_rxns_for_reactant(sp.id);
    if (rxn_classes == nullptr) {
      continue;
    }

    // go through all applicable reactants
    for (auto it: *rxn_classes) {
      const BNG::RxnClass* rxn_class = it.second;
      assert(rxn_class->reactants.size() == 2);

      species_id_t second_species_id;
      if (rxn_class->reactants[0] != sp.id) {
        second_species_id = rxn_class->reactants[0];
      }
      else {
        second_species_id = rxn_class->reactants[1];
      }
      const BNG::Species& sp2 = all_species.get(second_species_id);

      // we can use is_vol/is_surf for ALL_VOLUME_MOLECULES and ALL_SURFACE_MOLECULES
      if (sp.is_vol() || sp.id == all_molecules_species_id) {
        if (sp2.is_vol() || sp2.id == all_molecules_species_id) {
          has_vol_vol_rxn = true;
          sp.set_flag(BNG::SPECIES_FLAG_CAN_VOLVOL);
        }
        if (sp2.is_surf() || sp2.id == all_molecules_species_id) {
          sp.set_flag(BNG::SPECIES_FLAG_CAN_VOLSURF);
        }
        if (sp2.is_reactive_surface()) {
          sp.set_flag(BNG::SPECIES_FLAG_CAN_VOLWALL);
          if (sp.id == all_molecules_species_id || sp.id == all_volume_molecules_species_id) {
            // superclasses are processed first, this setting is then passed-on on all subsequent vol species
            all_vol_mols_can_react_with_surface = true;
          }
        }
      }

      if ((sp.is_surf() || sp.id == all_molecules_species_id)) {

        if (sp2.is_surf() || sp2.id == all_molecules_species_id) {
          sp.set_flag(BNG::SPECIES_FLAG_CAN_SURFSURF);
        }

        if (sp2.is_reactive_surface()) {
          sp.set_flag(BNG::SPECIES_FLAG_CAN_REGION_BORDER);
          if (sp.id == all_surface_molecules_species_id || sp.id == all_molecules_species_id) {
            all_surf_mols_can_react_with_surface = true;
          }
        }
      }
    }
  }

  // for mcell3 compatibility
  /* If there are no 3D molecules-reactants in the simulation
     set up the "use_expanded_list" flag to zero. */
  if (!has_vol_vol_rxn) {
    world->config.use_expanded_list = 0;
  }
}


MCell::wall_index_t MCell4Converter::convert_wall_and_add_to_geom_object(
    const API::GeometryObject& src_obj, const uint side,
    MCell::Partition& p, MCell::GeometryObject& dst_obj) {

  assert(src_obj.element_connections[side].size() == 3);

  // TODO LATER: there is really no reason to add walls in two steps,
  // can be simplified
  MCell::Wall& wall = p.add_uninitialized_wall(world->get_next_wall_id());

  wall.object_id = dst_obj.id;
  wall.object_index = dst_obj.index;
  wall.side = side;

  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    wall.vertex_indices[i] = src_obj.vertex_indices[src_obj.element_connections[side][i]];
  }

  wall.precompute_wall_constants(p);

  // add wall to subpartitions
  p.finalize_wall_creation(wall.index);
  dst_obj.wall_indices.push_back(wall.index);

  return wall.index;
}


MCell::region_index_t MCell4Converter::convert_surface_region(
    MCell::Partition& p,
    API::SurfaceRegion& surface_region, API::GeometryObject& o,
    MCell::GeometryObject& obj) {

  MCell::Region reg;
  reg.name = obj.name + "," + surface_region.name;
  reg.geometry_object_id = obj.id;

  // simply add all walls
  for (const int wall_in_object: surface_region.wall_indices) {
    wall_index_t wi = o.wall_indices[wall_in_object];
    reg.add_wall_to_walls_and_edges(wi, false);
  }

  if (is_set(surface_region.surface_class)) {
    assert(surface_region.surface_class->species_id != SPECIES_ID_INVALID);
    reg.species_id = surface_region.surface_class->species_id;
  }

  reg.init_surface_region_edges(p);

  reg.id = world->get_next_region_id();
  surface_region.region_id = reg.id;
  region_index_t ri = p.add_region_and_set_its_index(reg);
  return ri;
}


void MCell4Converter::convert_geometry_objects() {
  MCell::Partition& p = world->get_partition(0); // only partition 0 is supported for now

  for (std::shared_ptr<API::GeometryObject>& o: model->geometry_objects) {

    MCell::GeometryObject& obj = p.add_uninitialized_geometry_object(world->get_next_geometry_object_id());

    obj.name = o->name;
    obj.parent_name = ""; // empty, we do not really care
    o->geometry_object_id = obj.id;

    // vertices
    for (auto& v: o->vertex_list) {
      // add to partition and remember its index
      // must use rcp_length_unit to be identical to mcell4 with mdl
      vertex_index_t vi = p.add_or_find_geometry_vertex(Vec3(v[0], v[1], v[2]) * Vec3(world->config.rcp_length_unit));
      o->vertex_indices.push_back(vi);
    }

    // walls (validity of indices is checked in API::GeometryObject::check_semantics)
    for (size_t i = 0; i < o->element_connections.size(); i++) {
      wall_index_t wi = convert_wall_and_add_to_geom_object(*o, i, p, obj);
      o->wall_indices.push_back(wi);
    }

    // initialize edges
    obj.initialize_neighboring_walls_and_their_edges(p);

    vector<region_index_t> region_indices;
    // region "ALL"
    MCell::Region reg_all;
    reg_all.init_from_whole_geom_object(obj);
    reg_all.id = world->get_next_region_id();
    if (is_set(o->surface_class)) {
      assert(o->surface_class->species_id != SPECIES_ID_INVALID);
      reg_all.species_id = o->surface_class->species_id;
    }


    region_index_t ri_all = p.add_region_and_set_its_index(reg_all);
    region_indices.push_back(ri_all);

    // we must remember that this region belongs to our object
    obj.encompassing_region_id = ri_all;

    // regions from surface areas
    // mcell3 stores the regions in reverse, so we can too...
    // (maybe change this in the mcell3 converter)
    for (int k = (int)o->surface_regions.size() - 1; k >= 0; k--) {
      std::shared_ptr<SurfaceRegion> surface_region = o->surface_regions[k];
      region_index_t ri = convert_surface_region(p, *surface_region, *o, obj);
      region_indices.push_back(ri);
    }

    // also set regions for walls
    for (wall_index_t wi: obj.wall_indices) {
      Wall& w = p.get_wall(wi);

      for (region_index_t ri: region_indices) {
        MCell::Region& reg = p.get_region(ri);
        if (reg.walls_and_edges.count(wi) == 1) {
          w.regions.insert(ri);
        }
      }
    }

    // and assign surface class species to regions
  }
}


static MCell::RegionExprOperator convert_region_node_type(API::RegionNodeType t) {
  switch (t) {
    case API::RegionNodeType::LEAF_GEOMETRY_OBJECT:
    case API::RegionNodeType::LEAF_SURFACE_REGION:
      return MCell::RegionExprOperator::Leaf;

    case API::RegionNodeType::UNION:
      return MCell::RegionExprOperator::Union;

    case API::RegionNodeType::DIFFERENCE:
      return MCell::RegionExprOperator::Difference;

    case API::RegionNodeType::INTERSECT:
      return MCell::RegionExprOperator::Intersect;

    default:
      assert(false);
      return MCell::RegionExprOperator::Invalid;
  }
}


RegionExprNode* convert_region_expr_recursively(
    const shared_ptr<API::Region>& region,
    MCell::ReleaseEvent* rel_event
) {
  assert(is_set(region));
  if (region->node_type == RegionNodeType::LEAF_GEOMETRY_OBJECT) {
    shared_ptr<API::GeometryObject> geometry_object = dynamic_pointer_cast<API::GeometryObject>(region);
    assert(is_set(geometry_object));
    return rel_event->create_new_region_expr_node_leaf(geometry_object->name + MCell::REGION_ALL_SUFFIX_W_COMMA);
  }
  else if (region->node_type == RegionNodeType::LEAF_SURFACE_REGION) {
    shared_ptr<API::SurfaceRegion> surface_region = dynamic_pointer_cast<API::SurfaceRegion>(region);
    assert(is_set(surface_region));
    return rel_event->create_new_region_expr_node_leaf(surface_region->parent->name + "," + surface_region->name);
  }
  else {
    return rel_event->create_new_region_expr_node_op(
        convert_region_node_type(region->node_type),
        convert_region_expr_recursively(region->left_node, rel_event),
        convert_region_expr_recursively(region->right_node, rel_event)
    );
  }
}


void MCell4Converter::convert_region_expr(API::ReleaseSite& rel_site, MCell::ReleaseEvent* rel_event) {

  if (!is_set(rel_site.region)) {
    throw RuntimeError("Region for release site " + rel_site.name + " was not set.");
  }

  rel_event->region_expr_root = convert_region_expr_recursively(rel_site.region, rel_event);

  // also set llf and urb
  bool ok = Geometry::get_region_expr_bounding_box(world, rel_event->region_expr_root, rel_event->region_llf, rel_event->region_urb);
  if (!ok) {
    throw RuntimeError("Region for releases specified by " + rel_event->region_expr_root->region_name + " is not closed.");
  }
}


void MCell4Converter::convert_molecule_list(
    const std::vector<std::shared_ptr<MoleculeReleaseInfo>>& molecule_list,
    MCell::ReleaseEvent* rel_event) {

  for (auto& item: molecule_list) {
    MCell::SingleMoleculeReleaseInfo info;
    info.species_id = item->species->species_id;
    assert(item->location.size() == 3);
    info.pos.x = item->location[0] * world->config.rcp_length_unit;
    info.pos.y = item->location[1] * world->config.rcp_length_unit;
    info.pos.z = item->location[2] * world->config.rcp_length_unit;
    info.orientation = convert_orientation(item->orientation); // not set means random orientation for surf mols

    if (world->get_all_species().get(info.species_id).is_vol() &&
        rel_event->orientation != ORIENTATION_NONE && rel_event->orientation != ORIENTATION_NOT_SET) {
      throw ValueError(
          S(NAME_CLASS_RELEASE_SITE) + " releases a volume molecule but orientation is set.");
    }

    rel_event->molecule_list.push_back(info);
  }
}


void MCell4Converter::convert_release_events() {
  // only initial support without any geometries

  for (std::shared_ptr<API::ReleaseSite>& r: model->release_sites) {
    MCell::ReleaseEvent* rel_event = new ReleaseEvent(world);

    rel_event->release_site_name = r->name;

    if (!is_set(r->molecule_list)) {
      rel_event->species_id = r->species->species_id;
      rel_event->orientation = convert_orientation(r->orientation);

      // FIXME: this should belong in the ReleaseSite ctor
      if (world->get_all_species().get(rel_event->species_id).is_surf() &&
          rel_event->orientation != ORIENTATION_UP && rel_event->orientation != ORIENTATION_DOWN) {
        throw ValueError(
            S(NAME_CLASS_RELEASE_SITE) + " " + r->name +
            " releases a surface molecule but orientation is not set.");
      }

      if (world->get_all_species().get(rel_event->species_id).is_vol() &&
          rel_event->orientation != ORIENTATION_NONE && rel_event->orientation != ORIENTATION_NOT_SET) {
        throw ValueError(
            S(NAME_CLASS_RELEASE_SITE) + " " + r->name +
            " releases a volume molecule but orientation is set.");
      }
    }

    rel_event->delay = r->release_time / world->config.time_unit;

    // release pattern
    if (is_set(r->release_pattern)) {
      API::ReleasePattern rp = *r->release_pattern;
      rel_event->number_of_trains = rp.number_of_trains;
      rel_event->train_interval = rp.train_interval / world->config.time_unit;
      rel_event->train_duration = rp.train_duration / world->config.time_unit;
      rel_event->release_interval = rp.release_interval / world->config.time_unit;
    }

    // release_number_method
    if (is_set(r->number_to_release)) {
      rel_event->release_number_method = ReleaseNumberMethod::ConstNum;
      rel_event->release_number = r->number_to_release;
    }
    else if (is_set(r->density)) {
      rel_event->release_number_method = ReleaseNumberMethod::DensityNum;
      rel_event->concentration = r->density;
    }
    else if (is_set(r->molecule_list)) {
      rel_event->release_number_method = ReleaseNumberMethod::ConstNum;
      convert_molecule_list(r->molecule_list, rel_event);
    }
    else {
      throw RuntimeError(
          S("The only supported release number type now is constant number specified with ") + NAME_NUMBER_TO_RELEASE + "."
      );
    }

    // release_shape
    switch (r->shape) {
      case Shape::SPHERICAL:
        rel_event->release_shape = ReleaseShape::SPHERICAL;
        rel_event->location = r->location * world->config.rcp_length_unit;
        rel_event->diameter = r->site_diameter * world->config.rcp_length_unit;
        break;
      case Shape::REGION_EXPR: {
          rel_event->release_shape = ReleaseShape::REGION;
          convert_region_expr(*r, rel_event);
          bool ok = rel_event->initialize_walls_for_release();
          if (!ok) {
            throw RuntimeError("Only simple surface regions are supported now, error for " + r->name + ".");
          }
        }
        break;
      case Shape::LIST:
        rel_event->release_shape = ReleaseShape::LIST;
        rel_event->diameter = r->site_diameter * world->config.rcp_length_unit;
        break;
      default:
        // should be caught earlier
        throw RuntimeError("The only supported shape now is Spherical.");
    }

    rel_event->update_event_time_for_next_scheduled_time();
    world->scheduler.schedule_event(rel_event);
  }
}


MCell::MolOrRxnCountTerm MCell4Converter::convert_count_term_leaf_and_init_counting_flags(
    const std::shared_ptr<API::CountTerm> ct,
    const int sign
) {
  MCell::MolOrRxnCountTerm res;
  res.sign_in_expression = sign;

  assert(is_set(ct));
  assert(ct->node_type == API::ExprNodeType::LEAF);

  // where
  bool is_obj_not_surf_reg = false; // to silence compiler warning
  geometry_object_id_t obj_id = GEOMETRY_OBJECT_ID_INVALID;
  region_id_t reg_id = REGION_ID_INVALID;
  if (is_set(ct->region)) {
    if (ct->region->node_type == API::RegionNodeType::LEAF_GEOMETRY_OBJECT) {
      is_obj_not_surf_reg = true;
      obj_id = dynamic_pointer_cast<GeometryObject>(ct->region)->geometry_object_id;
    }
    else if (ct->region->node_type == API::RegionNodeType::LEAF_SURFACE_REGION) {
      is_obj_not_surf_reg = false;
      reg_id = dynamic_pointer_cast<SurfaceRegion>(ct->region)->region_id;
    }
    else {
      // already checked in check_semantics
      assert(false && "Invalid region node type.");
    }
  }

  // what
  if (is_set(ct->species)) {

    res.species_id = ct->species->species_id;
    BNG::Species& sp = world->get_all_species().get(res.species_id);

    res.orientation = convert_orientation(ct->orientation);

    if (is_set(ct->region)) {
      if (sp.is_vol()) {
        if (!is_obj_not_surf_reg) {
          throw ValueError("Counting surface molecules " + sp.name + " on a surface is not allowed.");
        }

        res.type = MCell::CountType::EnclosedInObject;
        res.geometry_object_id = obj_id;

        // set species flag
        sp.set_flag(BNG::SPECIES_FLAG_COUNT_ENCLOSED);
        sp.set_flag(BNG::SPECIES_FLAG_COUNT_CONTENTS);

        // and also mark the object that we are counting molecules inside
        world->get_geometry_object(res.geometry_object_id).is_counted_volume = true;
      }
      else {
        // surf mols
        if (is_obj_not_surf_reg) {
          // need to get the region of this object
          MCell::GeometryObject& obj = world->get_geometry_object(obj_id);
          assert(obj.encompassing_region_id != MCell::REGION_ID_INVALID);
          reg_id = obj.encompassing_region_id;
        }

        res.type = MCell::CountType::PresentOnSurfaceRegion;
        res.region_id = reg_id;
        sp.set_flag(BNG::SPECIES_FLAG_COUNT_CONTENTS);

        // these are only surface regions and there is no need to set that they are counted
      }
    }
    else {
      res.type = MCell::CountType::EnclosedInWorld;
      sp.set_flag(BNG::SPECIES_FLAG_COUNT_ENCLOSED);
    }
  }
  else if (is_set(ct->reaction_rule))
  {
    assert(!is_set(ct->reaction_rule->rev_rate));
    res.rxn_rule_id = ct->reaction_rule->fwd_rxn_rule_id;

    // is this a surface rxn? -> at least one of the reactants is a surface mol
    BNG::RxnRule* rxn = world->get_all_rxns().get_rxn_rule(res.rxn_rule_id);

    // need to set flag
    world->get_all_rxns().get_rxn_rule(res.rxn_rule_id)->set_is_counted();

    if (is_set(ct->region)) {
      if (!rxn->is_surf_rxn()) {
        if (is_obj_not_surf_reg) {
          res.type = MCell::CountType::RxnCountInObject;
          res.geometry_object_id = obj_id;

          world->get_geometry_object(res.geometry_object_id).is_counted_volume = true;
        }
      }
      else {
        if (is_obj_not_surf_reg) {
          // need to get the region of this object
          MCell::GeometryObject& obj = world->get_geometry_object(obj_id);
          assert(obj.encompassing_region_id != MCell::REGION_ID_INVALID);
          reg_id = obj.encompassing_region_id;
        }

        res.type = MCell::CountType::RxnCountOnSurfaceRegion;
        res.region_id = reg_id;

        // these are only surface regions and there is no need to set that they are counted
      }
    }
    else {
      res.type = MCell::CountType::RxnCountInWorld;
    }
  }
  else {
    assert(false);
  }

  return res;
}


void MCell4Converter::convert_count_terms_recursively(
    const std::shared_ptr<API::CountTerm> ct,
    const int sign,
    MCell::MolOrRxnCountInfo& info
) {
  assert(is_set(ct));

  if (ct->node_type == API::ExprNodeType::LEAF) {
    MCell::MolOrRxnCountTerm term = convert_count_term_leaf_and_init_counting_flags(ct, sign);
    info.terms.push_back(term);
  }
  else if (ct->node_type == API::ExprNodeType::ADD || ct->node_type == API::ExprNodeType::SUB) {
    convert_count_terms_recursively(ct->left_node, sign, info);

    int next_sign;
    if (ct->node_type == API::ExprNodeType::SUB) {
      next_sign = -sign;
    }
    else {
      next_sign = sign;
    }

    convert_count_terms_recursively(ct->right_node, next_sign, info);
  }
  else {
    // cannot really happen
    throw RuntimeError("Invalid node_type in CountTerm.");
  }
}


void MCell4Converter::convert_mol_or_rxn_count_events_and_init_counting_flags() {
  for (const std::shared_ptr<API::Count>& c: model->counts) {
    MCell::MolOrRxnCountEvent* count_event = new MCell::MolOrRxnCountEvent(world);

    count_event->event_time = 0;
    count_event->periodicity_interval = c->every_n_timesteps;

    // create buffer
    count_buffer_id_t buffer_id =
        world->create_count_buffer(c->filename, API::DEFAULT_COUNT_BUFFER_SIZE);

    MCell::MolOrRxnCountInfo info(buffer_id);

    // process count terms
    if (is_set(c->count_expression)) {
      convert_count_terms_recursively(c->count_expression, +1, info);
    }
    else {
      convert_count_terms_recursively(dynamic_pointer_cast<API::CountTerm>(c), +1, info);
    }

    // having multiple MolOrRxnCountInfo per MolOrRxnCountEvent
    // was useful for MCell3 conversion, however for pymcell4 each count is a separate event
    count_event->add_mol_count_info(info);
    world->scheduler.schedule_event(count_event);
  }
}


void MCell4Converter::convert_viz_output_events() {
  for (std::shared_ptr<API::VizOutput>& v: model->viz_outputs) {
    MCell::VizOutputEvent* viz_event = new VizOutputEvent(world);

    viz_event->event_time = 0.0;
    viz_event->periodicity_interval = v->every_n_timesteps;
    viz_event->viz_mode = convert_viz_mode(v->mode);
    viz_event->file_prefix_name = v->filename_prefix;

    for (std::shared_ptr<API::Species>& s: v->species_list) {
      viz_event->species_ids_to_visualize.insert(s->species_id);
    }

    world->scheduler.schedule_event(viz_event);
  }
}


void MCell4Converter::add_ctrl_c_termination_event() {
  MCell::PeriodicCallEvent* event = new PeriodicCallEvent(world);
  event->event_time = 0;
  event->periodicity_interval = 1;
  event->function_ptr = check_ctrl_c;

  world->scheduler.schedule_event(event);
}

} // namespace API
} // namespace MCell
