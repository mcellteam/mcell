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
#include "rng.h"
#include "bng/bng.h"

#include "api/mcell.h"

using namespace std;

namespace MCell {
namespace API {

static orientation_t convert_orientation(const Orientation o) {
  switch (o) {
    case Orientation::Down:
      return ORIENTATION_DOWN;
    case Orientation::None:
      return ORIENTATION_NONE;
    case Orientation::Up:
      return ORIENTATION_UP;
    case Orientation::NotSet:
      return ORIENTATION_NOT_SET;
    default:
      throw ValueError("Invalid Orientation value " + to_string((int)o) + ".");
  }
}


static viz_mode_t convert_viz_mode(const VizMode m) {
  switch (m) {
    case VizMode::Ascii:
      return ASCII_MODE;
    case VizMode::Cellblender:
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

  convert_rxns();

  // at this point, we need to create the first (and for now the only) partition
  // create initial partition with center at 0,0,0
  partition_id_t index = world->add_partition(Vec3(0, 0, 0));
  assert(index == PARTITION_ID_INITIAL);

  convert_geometry_objects();

  convert_release_events();

  convert_viz_output_events();
}


void MCell4Converter::convert_simulation_setup() {
  const API::Config& config = model->config;

  if (config.microscopic_reversibility) {
    throw ValueError("Setting config.microscopic_reversibility to True is not supported yet.");
  }

  world->total_iterations = config.total_iterations_hint;
  world->config.time_unit = config.time_step;

  float_t grid_density = config.surface_grid_density;
  world->config.grid_density = grid_density;

  float_t length_unit = 1/sqrt_f(config.surface_grid_density);
  world->config.length_unit = length_unit;

  if (is_set(config.interaction_radius)) {
    world->config.rx_radius_3d = config.interaction_radius / length_unit;
  }
  else {
    world->config.rx_radius_3d = (1.0 / sqrt_f(MY_PI * grid_density)) / length_unit;
  }

  // TODO CHECK: converting to internal length units unlike as in mcell3
  world->config.vacancy_search_dist2 = config.vacancy_search_distance / length_unit;

  world->config.randomize_smol_pos = !config.center_molecules_on_grid;

  world->seed_seq = config.seed;
  rng_init(&world->rng, world->seed_seq);

  world->config.partition_edge_length = config.partition_dimension / length_unit;
  int num_subparts = config.partition_dimension / config.subpartition_dimension;
  assert(num_subparts > 0);
  if (num_subparts % 2 == 1) {
    // the number of subparts must be even
    num_subparts++;
  }
  world->config.subpartitions_per_partition_dimension = num_subparts;

  // this option in MCell3 was removed in MCell4
  world->config.use_expanded_list = true;

  // compute other constants
  world->config.init();
}


void MCell4Converter::convert_species() {
  for (std::shared_ptr<API::Species>& s: model->species) {
    BNG::Species new_species;
    new_species.name = s->name;

    bool is_vol;
    if (is_set(s->diffusion_constant_3d)) {
      new_species.D = s->diffusion_constant_3d;
      new_species.set_is_vol();
      is_vol = true;
    }
    else if (is_set(s->diffusion_constant_2d)) {
      new_species.D = s->diffusion_constant_2d;
      new_species.set_is_surf();
      is_vol = true;
    }
    else {
      throw ValueError("Neither diffusion_constant_2d nor diffusion_constant_3d was set.");
    }
	
		// TODO: set all flags that are used in MCell4

    new_species.update_space_and_time_step(world->config.time_unit, world->config.length_unit);

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
  }
}


BNG::ComponentType MCell4Converter::convert_component_type(API::ComponentType& ct) {
  throw RuntimeError("Components are not supported yet");
}


BNG::ComponentInstance MCell4Converter::convert_component_instance(API::ComponentInstance& ci) {
  throw RuntimeError("Components are not supported yet");
}


BNG::MolType MCell4Converter::convert_molecule_type(API::MoleculeType& mt) {
  BNG::MolType res;

  res.name = mt.name;
  if (!mt.components.empty()) {
    throw RuntimeError("Components are not supported yet");
  }

  return res;
}


BNG::MolInstance MCell4Converter::convert_molecule_instance(API::MoleculeInstance& mi) {
  BNG::MolInstance res;

  BNG::MolType mt = convert_molecule_type(*mi.molecule_type);
  res.mol_type_id = world->bng_engine.get_data().find_or_add_molecule_type(mt);

  if (!mi.components.empty()) {
    throw RuntimeError("Components are not supported yet");
  }

  return res;
}


BNG::CplxInstance MCell4Converter::convert_complex_instance(API::ComplexInstance& inst) {
  BNG::CplxInstance res;

  for (std::shared_ptr<API::MoleculeInstance>& m: inst.molecule_instances) {
    BNG::MolInstance mi = convert_molecule_instance(*m);

    res.mol_instances.push_back(mi);
  }

  res.set_orientation(convert_orientation(inst.orientation));
  res.finalize();
  return res;
}


void MCell4Converter::convert_rxns() {
  for (std::shared_ptr<API::ReactionRule>& r: model->reaction_rules) {

    bool is_reversible = is_set(r->rev_rate);
    if (is_reversible && is_set(r->name)) {
        assert(false); // checked with semantic check
    }

    BNG::RxnRule rxn;

    if (!is_reversible && is_set(r->name)) {
      rxn.name = r->name;

    }
    rxn.type = BNG::RxnType::Standard;
    rxn.rate_constant = r->fwd_rate;

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

    // TODO: variable rates

    // reverse reaction
    BNG::RxnRule rxn_rev;
    if (is_reversible) {
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


wall_index_t MCell4Converter::convert_wall_and_add_to_geom_object(
    const API::GeometryObject& src_obj, const uint side,
    MCell::Partition& p, MCell::GeometryObject& dst_obj) {

  assert(src_obj.element_connections[side].size() == 3);

  // TODO LATER: there is really no reason to add walls in two steps,
  // can be simplified
  Wall& wall = p.add_uninitialized_wall(world->get_next_wall_id());

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


void MCell4Converter::convert_geometry_objects() {
  MCell::Partition& p = world->get_partition(0); // only partition 0 is supported for now

  for (std::shared_ptr<API::GeometryObject>& o: model->geometry_objects) {

    MCell::GeometryObject& obj = p.add_uninitialized_geometry_object(world->get_next_geometry_object_id());

    obj.name = o->name;
    obj.parent_name = ""; // empty, we do not really care

    // vertices
    for (auto& v: o->vertex_list) {
      // add to partition and remember its index
      vertex_index_t vi = p.add_or_find_geometry_vertex(Vec3(v[0], v[1], v[2]) / Vec3(world->config.length_unit));
      o->vertex_indices.push_back(vi);
    }

    // walls (validity of indices is checked in API::GeometryObject::check_semantics)
    for (size_t i = 0; i < o->element_connections.size(); i++) {
      wall_index_t wi = convert_wall_and_add_to_geom_object(*o, i, p, obj);
      o->wall_indices.push_back(wi);
    }

    // initialize edges
    obj.initialize_neighboring_walls_and_their_edges(p);

    // regions
    assert(o->surface_regions.empty() && "TODO");
  }
}


void MCell4Converter::convert_release_events() {
  // only initial support without any geometries

  for (std::shared_ptr<API::ReleaseSite>& r: model->release_sites) {
    MCell::ReleaseEvent* rel_event = new ReleaseEvent(world);

    // release patterns are not supported yet
    rel_event->release_site_name = r->name;
    rel_event->event_time = 0.0;
    rel_event->actual_release_time = 0.0;
    rel_event->species_id = r->species->species_id;
    rel_event->orientation = convert_orientation(r->initial_orientation);

    switch (r->shape) {
      case Shape::Spherical:
        rel_event->release_shape = ReleaseShape::SPHERICAL;
        rel_event->location = r->location / world->config.length_unit;
        rel_event->diameter = r->site_diameter;
        break;
      default:
        // should be caught earlier
        throw RuntimeError("The only supported shape now is Spherical.");
    }

    if (is_set(r->number_to_release)) {
      rel_event->release_number_method = ReleaseNumberMethod::ConstNum;
      rel_event->release_number = r->number_to_release;
    }
    else {
      throw RuntimeError(
          "The only supported release number type now is constant number specified with 'number_to_release'."
      );
    }

    world->scheduler.schedule_event(rel_event);
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


} // namespace API
} // namespace MCell