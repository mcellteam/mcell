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

  convert_elementary_molecule_types();

  convert_species();

  convert_surface_classes();

  convert_rxns();

  // find out whether we have a vol vol rxn for all current species (for mcell3 compatibility)
  // in BNG mode this finds a reaction as well
  // this is an optimization to tell that we don't need to check the
  // surroundings of subpartitions in case when there are no vol-vol rxns
  // WARNING: must be done before geometry objects are converted because wall_subparts_collision_test
  // depends on it
  if (!world->get_all_rxns().has_bimol_vol_rxns()) {
    // the default is true (or read from user)
    world->config.use_expanded_list = false;
  }

  // at this point, we need to create the first (and for now the only) partition
  // create initial partition with center at 0,0,0
  partition_id_t index = world->add_partition(world->config.partition0_llf);
  assert(index == PARTITION_ID_INITIAL);

  convert_geometry_objects();

  // uses random generator state
  Geometry::check_for_overlapped_walls(world);

  // we need to schedule the initial release for surfaces before the other releases
  world->create_initial_surface_region_release_event();

  convert_release_events();

  // diffusion events must be created after reactions and release events
  // since they may define the initial species in when BNGL reactions are used
  world->create_diffusion_events();

  convert_mol_or_rxn_count_events_and_init_counting_flags();

  convert_viz_output_events();

  // must be called after release events define species because they don't have to be
  // defined, also after count events were defined
  init_rxn_related_flags();

  add_ctrl_c_termination_event();

  // some general checks
  if (world->config.rx_radius_3d >= world->config.subpartition_edge_length) {
    throw ValueError(S("Reaction radius ") + to_string(world->config.rx_radius_3d * world->config.length_unit) +
        " is larger than subpartition edge length " + to_string(world->config.subpartition_edge_length * world->config.length_unit) + ".");
  }
}


void MCell4Converter::get_geometry_bounding_box(Vec3& llf, Vec3& urb) {

  llf = Vec3(FLT_GIGANTIC);
  urb = Vec3(-FLT_GIGANTIC);

  for (std::shared_ptr<API::GeometryObject>& o: model->geometry_objects) {
    // go through all vertices
    for (auto& vert: o->vertex_list) {
      Vec3 v(vert);

      if (v.x < llf.x) {
        llf.x = v.x;
      }
      if (v.x > urb.x) {
        urb.x = v.x;
      }
      if (v.y < llf.y) {
        llf.y = v.y;
      }
      if (v.y > urb.y) {
        urb.y = v.y;
      }
      if (v.z < llf.z) {
        llf.z = v.z;
      }
      if (v.z > urb.z) {
        urb.z = v.z;
      }

    }
  }
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

  world->config.initial_seed = config.seed;
  rng_init(&world->rng, world->config.initial_seed);


  Vec3 llf, urb;
  get_geometry_bounding_box(llf, urb);
  Vec3 llf_w_margin = llf - Vec3(PARTITION_EDGE_EXTRA_MARGIN_UM);
  Vec3 urb_w_margin = urb + Vec3(PARTITION_EDGE_EXTRA_MARGIN_UM);

  float_t llf_partition_dimension_diff = 0;

  // TODO: should we be modifying the values user provided?

  float_t auto_partition_dimension = max3(urb_w_margin) - min3(llf_w_margin);
  if (!is_set(config.initial_partition_origin) && auto_partition_dimension > config.partition_dimension) {
    cout <<
        "Value of " << NAME_CLASS_MODEL << "." << NAME_CONFIG << "." << NAME_PARTITION_DIMENSION <<
        " " << config.partition_dimension <<
        " is smaller than the automatically determined value " << auto_partition_dimension <<
        ", using the automatic value.\n";
    world->config.partition_edge_length = auto_partition_dimension / length_unit;
    // we need to move the origin, if specified, by half of this increment
    llf_partition_dimension_diff =
        -(world->config.partition_edge_length - (config.partition_dimension / length_unit)) / 2;
  }
  else if (is_set(config.initial_partition_origin) &&
      (config.initial_partition_origin[0] > llf_w_margin.x ||
       config.initial_partition_origin[1] > llf_w_margin.y ||
       config.initial_partition_origin[2] > llf_w_margin.z ||
       config.initial_partition_origin[0] + config.partition_dimension < urb_w_margin.x ||
       config.initial_partition_origin[1] + config.partition_dimension < urb_w_margin.y ||
       config.initial_partition_origin[2] + config.partition_dimension < urb_w_margin.z
      )
  ) {
    Vec3 origin(config.initial_partition_origin);
    Vec3 llf_diff(0);
    if (origin.x > llf_w_margin.x) {
      llf_diff.x = origin.x - llf_w_margin.x;
    }
    if (origin.y > llf_w_margin.y) {
      llf_diff.y = origin.y - llf_w_margin.y;
    }
    if (origin.z > llf_w_margin.z) {
      llf_diff.z = origin.z - llf_w_margin.z;
    }
    llf_partition_dimension_diff = max3(llf_diff);

    Vec3 opposite(Vec3(config.initial_partition_origin) + Vec3(config.partition_dimension));
    Vec3 urb_diff(0);
    if (opposite.x < urb_w_margin.x) {
      urb_diff.x = urb_w_margin.x - opposite.x;
    }
    if (opposite.y < urb_w_margin.y) {
      urb_diff.y = urb_w_margin.y - opposite.y;
    }
    if (opposite.z < urb_w_margin.z) {
      urb_diff.z = urb_w_margin.z - opposite.z;
    }
    float_t max_urb_diff = max3(urb_diff);

    world->config.partition_edge_length =
        (config.partition_dimension + llf_partition_dimension_diff + max_urb_diff)/ length_unit;
    // we need to move the origin, if specified, by half of this increment

    Vec3 new_origin = origin - Vec3(llf_partition_dimension_diff);

    cout <<
        "Value of " << NAME_INITIAL_PARTITION_ORIGIN << " " << origin <<
        " does not provide enough margin"
        " for model's geometry bounding box lower, left, front point " << llf << " "
        " and upper, right, back " << urb << "."
        " Moving the origin to " << new_origin << " and increasing " <<
        NAME_PARTITION_DIMENSION << " to " <<
        world->config.partition_edge_length * length_unit << ".\n";
  }
  else {
    world->config.partition_edge_length = config.partition_dimension / length_unit;
  }

  if (is_set(config.initial_partition_origin)) {
    // origin set manually
    world->config.partition0_llf =
        (Vec3(config.initial_partition_origin) - Vec3(llf_partition_dimension_diff))/ Vec3(length_unit);
  }
  else {
    // place the partition to the center
    world->config.partition0_llf = -Vec3(world->config.partition_edge_length) / Vec3(2);
  }

  // align the origin to a multiple of subpartition length
  float_t sp_len = config.subpartition_dimension / length_unit;

  uint tentative_subparts = world->config.partition_edge_length / sp_len;
  if (tentative_subparts > MAX_SUBPARTS_PER_PARTITION) {
    cout <<
      "Approximate number of subpartitions " << tentative_subparts <<
      " is too high, lowering it to a limit of " << MAX_SUBPARTS_PER_PARTITION << ".\n";
    sp_len = world->config.partition_edge_length / MAX_SUBPARTS_PER_PARTITION;
  }

  Vec3 orig_origin = world->config.partition0_llf;
  world->config.partition0_llf =
      floor_to_multiple_allow_negative(orig_origin, sp_len);

  // enlarge the partition by size we moved it in order to be aligned
  Vec3 llf_moved = orig_origin - world->config.partition0_llf;
  float_t partition_edge_length_enlarged = world->config.partition_edge_length + max3(llf_moved);
  world->config.partition_edge_length = ceil_to_multiple(partition_edge_length_enlarged, sp_len);

  world->config.num_subpartitions_per_partition =
      round_f(world->config.partition_edge_length / sp_len);


  // this option in MCell3 was removed in MCell4
  world->config.use_expanded_list = true;

  // compute other constants
  world->config.init();
}


BNG::mol_type_id_t MCell4Converter::convert_elementary_molecule_type(API::ElementaryMoleculeType& api_mt, const bool in_rxn) {
  if (api_mt.mol_type_id != BNG::MOL_TYPE_ID_INVALID) {
    // already converted
    return api_mt.mol_type_id;
  }


  BNG::BNGData& bng_data = world->bng_engine.get_data();

  BNG::MolType bng_mt;

  bng_mt.name = api_mt.name;

  if (!in_rxn) {
    if (is_set(api_mt.diffusion_constant_2d)) {
      bng_mt.set_is_surf();
      bng_mt.D = api_mt.diffusion_constant_2d;
    }
    else if (is_set(api_mt.diffusion_constant_3d)) {
      bng_mt.set_is_vol();
      bng_mt.D = api_mt.diffusion_constant_3d;
    }
    else {
      throw RuntimeError(S("Diffusion constant for ") + NAME_CLASS_ELEMENTARY_MOLECULE_TYPE + " was not set.");
    }

    if (is_set(api_mt.custom_time_step)) {
      bng_mt.custom_time_step = api_mt.custom_time_step;
    }
    else if (is_set(api_mt.custom_space_step)) {
      bng_mt.custom_space_step = api_mt.custom_space_step;
    }
  }

  // components
  for (std::shared_ptr<API::ComponentType> api_ct: api_mt.components) {
    BNG::ComponentType bng_ct;

    bng_ct.name = api_ct->name;

    for (const string& state: api_ct->states) {
      bng_ct.allowed_state_ids.insert_unique(bng_data.find_or_add_state_name(state));
    }

    bng_mt.component_type_ids.push_back(bng_data.find_or_add_component_type(bng_ct));
  }

  return bng_data.find_or_add_molecule_type(bng_mt);
}


void MCell4Converter::convert_elementary_molecule_types() {
  BNG::BNGData& bng_data = world->bng_engine.get_data();

  for (std::shared_ptr<API::ElementaryMoleculeType>& api_mt: model->elementary_molecule_types) {
    convert_elementary_molecule_type(*api_mt);
  }
}


void MCell4Converter::convert_species() {
  for (std::shared_ptr<API::Species>& s: model->species) {
    BNG::Species new_species(world->bng_engine.get_data());
    new_species.name = s->name;

    if (is_set(s->custom_time_step)) {
      new_species.custom_time_step = s->custom_time_step;
    }
    else if (is_set(s->custom_space_step)) {
      new_species.custom_space_step = s->custom_space_step;
    }

    bool is_vol = false;
    if (is_set(s->diffusion_constant_3d)) {
      is_vol = true;
      new_species.D = s->diffusion_constant_3d;
      new_species.set_is_vol();
      new_species.update_space_and_time_step(world->bng_engine.get_config());
    }
    else if (is_set(s->diffusion_constant_2d)) {
      is_vol = false;
      new_species.D = s->diffusion_constant_2d;
      new_species.set_is_surf();
      new_species.update_space_and_time_step(world->bng_engine.get_config());
    }
    else if (is_species_superclass(new_species.name)) {
      is_vol = new_species.name != ALL_SURFACE_MOLECULES;
      // these values are not really used, they are initialized for comparisons
      new_species.D = 0;
      new_species.space_step = 0;
      new_species.time_step = 0;
    }
    else {
      throw ValueError("Neither diffusion_constant_2d nor diffusion_constant_3d was set.");
    }
	
    new_species.set_flag(BNG::SPECIES_FLAG_CANT_INITIATE, s->target_only); // default is false

    // we must add a complex instance as the single molecule type in the new species
    // define a molecule type with no components
    BNG::MolType mol_type;
    mol_type.name = new_species.name; // name of the mol type is the same as for our species
    mol_type.D = new_species.D; // we must also set the diffusion constant - simply inherit from this simple species
    if (is_vol) {
      mol_type.set_is_vol();
    }
    else {
      mol_type.set_is_surf();
    }
    BNG::mol_type_id_t mol_type_id = world->bng_engine.get_data().find_or_add_molecule_type(mol_type);

    BNG::MolInstance mol_inst;
    mol_inst.mol_type_id = mol_type_id;

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

  BNG::RxnRule rxn(&world->bng_engine.get_data());

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
  rxn.base_rate_constant = FLT_GIGANTIC;

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
    BNG::Species sc_species(world->bng_engine.get_data());
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
    mol_type.set_is_reactive_surface();
    BNG::mol_type_id_t mol_type_id = world->bng_engine.get_data().find_or_add_molecule_type(mol_type);

    BNG::MolInstance mol_inst;
    mol_inst.mol_type_id = mol_type_id;
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


BNG::component_type_id_t MCell4Converter::convert_component_type(API::ComponentType& api_ct) {
  // component types are identified by their name and set of allowed states, not just by their name
  BNG::ComponentType bng_ct;

  bng_ct.name = api_ct.name;

  for (string& s: api_ct.states) {
    bng_ct.allowed_state_ids.insert( world->bng_engine.get_data().find_or_add_state_name(s) );
  }

  return world->bng_engine.get_data().find_or_add_component_type(bng_ct);
}


BNG::ComponentInstance MCell4Converter::convert_component_instance(API::ComponentInstance& api_ci) {

  BNG::ComponentInstance res(convert_component_type(*api_ci.component_type));

  if (api_ci.state == STATE_UNSET) {
    res.state_id = BNG::STATE_ID_DONT_CARE;
  }
  else {
    res.state_id = world->bng_engine.get_data().find_state_id(api_ci.state);
    assert(res.state_id != BNG::STATE_ID_INVALID);
  }

  if (api_ci.bond == BOND_BOUND) {
    res.bond_value = BNG::BOND_VALUE_BOUND;
  }
  else if (api_ci.bond == BOND_UNBOUND) {
    res.bond_value = BNG::BOND_VALUE_UNBOUND;
  }
  else if (api_ci.bond == BOND_ANY) {
    res.bond_value = BNG::BOND_VALUE_ANY;
  }
  else {
    res.bond_value = api_ci.bond;
    assert(res.bond_value != BNG::BOND_VALUE_INVALID);
  }

  return res;
}


BNG::MolInstance MCell4Converter::convert_molecule_instance(API::ElementaryMoleculeInstance& mi, const bool in_rxn) {
  BNG::MolInstance res;

  res.mol_type_id = convert_elementary_molecule_type(*mi.elementary_molecule_type, in_rxn);

  for (std::shared_ptr<API::ComponentInstance>& api_ci: mi.components) {
    res.component_instances.push_back(convert_component_instance(*api_ci));
  }

  // we must also copy flags from the mol type
  res.finalize_flags_and_sort_components(world->bng_engine.get_data());

  return res;
}


BNG::CplxInstance MCell4Converter::convert_complex_instance(API::ComplexInstance& inst, const bool in_rxn) {
  // create a temporary cplx instance that we will use for search
  BNG::CplxInstance cplx_inst(&world->bng_engine.get_data());

  for (std::shared_ptr<API::ElementaryMoleculeInstance>& m: inst.elementary_molecule_instances) {
    BNG::MolInstance mi = convert_molecule_instance(*m, in_rxn);

    cplx_inst.mol_instances.push_back(mi);
  }
  orientation_t orient = convert_orientation(inst.orientation, true);
  cplx_inst.set_orientation(orient);
  cplx_inst.finalize();

  if (!in_rxn) {
    // we need to find or add existing species that we match
    species_id_t species_id = world->get_all_species().find_full_match(cplx_inst);

    if (species_id == SPECIES_ID_INVALID) {
      // FIXME: can we have just one method to compare the cplx and create
      // species if needed?
      // cplx inst should clearly identify the species...
      species_id = world->get_all_species().find_or_add(
        BNG::Species(cplx_inst, world->bng_engine.get_data(), world->bng_engine.get_config())
      );
    }
    assert(species_id != SPECIES_ID_INVALID);
    return world->bng_engine.create_cplx_instance_for_species(species_id, orient);
  }
  else {
    return cplx_inst;
  }
}


void MCell4Converter::convert_rxns() {
  for (std::shared_ptr<API::ReactionRule>& r: model->reaction_rules) {

    bool is_reversible = is_set(r->rev_rate);

    BNG::RxnRule rxn(&world->bng_engine.get_data());

    if (is_set(r->name)) {
      rxn.name = r->name;
    }
    rxn.type = BNG::RxnType::Standard;

    if (is_set(r->fwd_rate)) {
      rxn.base_rate_constant = r->fwd_rate;
    }
    else {
      assert(!r->variable_rate.empty());
      for (auto& time_and_rate: r->variable_rate) {
        assert(time_and_rate.size() == 2);
        BNG::RxnRateInfo info;
        info.time = time_and_rate[0] / world->config.time_unit;
        info.rate_constant = time_and_rate[1];
        rxn.base_variable_rates.push_back(info);
      }
      rxn.update_variable_rxn_rate(0, nullptr);
    }

    for (std::shared_ptr<API::ComplexInstance>& rinst: r->reactants) {
      // convert to BNG::ComplexInstance using existing or new BNG::molecule_id

      BNG::CplxInstance reactant = convert_complex_instance(*rinst, true);
      rxn.append_reactant(reactant);
    }

    for (std::shared_ptr<API::ComplexInstance>& pinst: r->products) {
      // convert to BNG::ComplexInstance using existing or new BNG::molecule_id

      BNG::CplxInstance product = convert_complex_instance(*pinst, true);
      rxn.append_product(product);
    }

    // reverse reaction
    BNG::RxnRule rxn_rev(&world->bng_engine.get_data());
    if (is_reversible) {
      rxn_rev.type = BNG::RxnType::Standard;
      if (is_set(r->rev_name)) {
        rxn.name = r->rev_name;
      }
      rxn_rev.base_rate_constant = r->rev_rate;
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


void MCell4Converter::init_rxn_related_flags() {
  BNG::SpeciesContainer& all_species = world->get_all_species();
  BNG::RxnContainer& all_rxns = world->get_all_rxns();

  species_id_t all_molecules_species_id = all_species.get_all_molecules_species_id();
  species_id_t all_volume_molecules_species_id = all_species.get_all_volume_molecules_species_id();
  species_id_t all_surface_molecules_species_id = all_species.get_all_surface_molecules_species_id();

  all_rxns.update_all_mols_flags();

  // needed to update counting flags
  world->get_all_species().recompute_species_flags(world->get_all_rxns());
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

void MCell4Converter::convert_initial_surface_releases(
    const std::vector<std::shared_ptr<API::InitialSurfaceRelease>>& api_releases,
    std::vector<MCell::InitialRegionMolecules>& mcell_releases
) {
  for (auto api_rel: api_releases) {
    species_id_t species_id = api_rel->species->species_id;
    orientation_t orientation = convert_orientation(api_rel->orientation);

    if (is_set(api_rel->number_to_release)) {
      mcell_releases.push_back(
          InitialRegionMolecules(species_id, orientation, true, (uint)api_rel->number_to_release)
      );
    }
    else if (is_set(api_rel->density)) {
      mcell_releases.push_back(
          InitialRegionMolecules(species_id, orientation, false, (float_t)api_rel->density)
      );
    }
    else {
      assert(false);
    }
  }
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
  if (is_set(surface_region.initial_surface_releases)) {
    convert_initial_surface_releases(
        surface_region.initial_surface_releases,
        reg.initial_region_molecules
    );
  }

  reg.init_surface_region_edges(p);

  reg.id = world->get_next_region_id();
  surface_region.region_id = reg.id;
  region_index_t ri = p.add_region_and_set_its_index(reg);
  return ri;
}


void MCell4Converter::convert_geometry_objects() {
  MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL); // only partition 0 is supported for now

  for (std::shared_ptr<API::GeometryObject>& o: model->geometry_objects) {

    MCell::GeometryObject& obj = p.add_uninitialized_geometry_object(world->get_next_geometry_object_id());

    obj.name = o->name;
    obj.parent_name = ""; // empty, we do not really care
    o->geometry_object_id = obj.id;

    // vertices
    // remember the "offset" from the first vertex in the target partition
    o->first_vertex_index = p.get_geometry_vertex_count();
    for (auto& v: o->vertex_list) {
      // add to partition and remember its index
      // must use rcp_length_unit to be identical to mcell4 with mdl
      vertex_index_t vi =
          p.add_geometry_vertex(Vec3(v[0], v[1], v[2]) * Vec3(world->config.rcp_length_unit));
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
    if (is_set(o->initial_surface_releases)) {
      convert_initial_surface_releases(o->initial_surface_releases, reg_all.initial_region_molecules);
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


RegionExprNode* MCell4Converter::convert_region_expr_recursively(
    const shared_ptr<API::Region>& region,
    MCell::ReleaseEvent* rel_event
) {
  assert(is_set(region));
  if (region->node_type == RegionNodeType::LEAF_GEOMETRY_OBJECT) {
    shared_ptr<API::GeometryObject> geometry_object = dynamic_pointer_cast<API::GeometryObject>(region);
    assert(is_set(geometry_object));
    const MCell::Region* reg = world->find_region_by_name(geometry_object->name + MCell::REGION_ALL_SUFFIX_W_COMMA);
    release_assert(reg != nullptr);
    return rel_event->create_new_region_expr_node_leaf(reg->id);
  }
  else if (region->node_type == RegionNodeType::LEAF_SURFACE_REGION) {
    shared_ptr<API::SurfaceRegion> surface_region = dynamic_pointer_cast<API::SurfaceRegion>(region);
    assert(is_set(surface_region));
    const MCell::Region* reg = world->find_region_by_name(surface_region->parent->name + "," + surface_region->name);
    release_assert(reg != nullptr);
    return rel_event->create_new_region_expr_node_leaf(reg->id);
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
  // TODO: this does not check anything yet
  bool ok = Geometry::compute_region_expr_bounding_box(world, rel_event->region_expr_root, rel_event->region_llf, rel_event->region_urb);
  if (!ok) {
    throw RuntimeError("Region for releases specified by " + rel_event->region_expr_root->to_string(world) + " is not closed.");
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

      if (is_set(r->species)) {
        rel_event->species_id = r->species->species_id;
      }
      else if (is_set(r->complex_instance)) {
        // we need to define species for our complex instance
        BNG::Species s = BNG::Species(
            convert_complex_instance(*r->complex_instance),
            world->bng_engine.get_data(),
            world->bng_engine.get_config()
        );
        rel_event->species_id = world->bng_engine.get_all_species().find_or_add(s);
      }
      else {
        assert(false);
      }

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
    else if (is_set(r->concentration)) {
      rel_event->release_number_method = ReleaseNumberMethod::ConcentrationNum;
      rel_event->concentration = r->concentration;
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
            throw RuntimeError("Only simple surface regions are supported for surface releases currently, error for " + r->name + ".");
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

        res.type = MCell::CountType::EnclosedInVolumeRegion;
        res.geometry_object_id = obj_id;

        // set species flag
        sp.set_flag(BNG::SPECIES_FLAG_COUNT_ENCLOSED);
        sp.set_flag(BNG::SPECIES_FLAG_COUNT_CONTENTS);

        sp.set_flag(BNG::SPECIES_FLAG_NEEDS_COUNTED_VOLUME);

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
    BNG::RxnRule* rxn = world->get_all_rxns().get(res.rxn_rule_id);

    if (is_set(ct->region)) {
      if (!rxn->is_surf_rxn()) {
        // volume reaction
        if (is_obj_not_surf_reg) {
          res.type = MCell::CountType::RxnCountInVolumeRegion;
          rxn->set_is_counted_in_volume_regions();

          res.geometry_object_id = obj_id;
          world->get_geometry_object(res.geometry_object_id).is_counted_volume = true;
        }
        else {
          throw RuntimeError("Cannot count volume reaction " + rxn->name + " on surface " +
              ct->region->name + ".");
        }
      }
      else {
        // surface reaction
        res.type = MCell::CountType::RxnCountOnSurfaceRegion;
        rxn->set_is_counted_on_surface_regions();

        if (is_obj_not_surf_reg) {
          // need to get the region of this object
          MCell::GeometryObject& obj = world->get_geometry_object(obj_id);
          assert(obj.encompassing_region_id != MCell::REGION_ID_INVALID);
          reg_id = obj.encompassing_region_id;
        }

        res.region_id = reg_id;

        // these are only surface regions and there is no need to set that they are counted
      }
    }
    else {
      res.type = MCell::CountType::RxnCountInWorld;
      rxn->set_is_counted_in_world();
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

    info.multiplier = c->multiplier;

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

    if (is_set(v->species_list)) {
      for (std::shared_ptr<API::Species>& s: v->species_list) {
        viz_event->species_ids_to_visualize.insert(s->species_id);
      }
    }
    else if (v->all_species) {
      viz_event->visualize_all_species = true;
    }
    else {
      // empty, files will be created but there will be nothing to visualize
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
