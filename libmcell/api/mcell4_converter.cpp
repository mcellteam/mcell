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
#include "rng.h"

#include "api/config.h"
#include "api/species.h"

namespace MCell {
namespace API {


void MCell4Converter::convert(Model* model_, World* world_) {
  model = model_;
  world = world_;

  convert_simulation_setup();
  convert_species();
}


void MCell4Converter::convert_simulation_setup() {
  const Config& config = model->config;

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
    world->config.rx_radius_3d = 1.0 / sqrt_f(MY_PI * grid_density);
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
    if (is_set(s->diffusion_constant_2d)) {
      new_species.D = s->diffusion_constant_2d;
      new_species.set_is_vol();
    }
    else if (is_set(s->diffusion_constant_3d)) {
      new_species.D = s->diffusion_constant_3d;
      new_species.set_is_surf();
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
    if (new_species.is_vol()) {
      mol_inst.set_is_vol();
    }
    else if (new_species.is_surf()) {
      mol_inst.set_is_surf();
    }
    else {
      assert(false);
    }

    new_species.mol_instances.push_back(mol_inst);

    new_species.finalize();
    species_id_t new_species_id = world->get_all_species().find_or_add(new_species);

    // remember which species we created
    s->species_id = new_species_id;
  }
}


} // namespace API
} // namespace MCell
