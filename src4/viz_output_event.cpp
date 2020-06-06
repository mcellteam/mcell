/******************************************************************************
 *
 * Copyright (C) 2019 by
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


#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <errno.h>

#include "logging.h"
#include "mem_util.h"
#include "util.h"
#include "mcell_structs.h"

#include "viz_output_event.h"
#include "geometry.h"
#include "world.h"
#include "datamodel_defines.h"

#include "geometry_utils.inc"

using namespace std;

namespace MCell {

void VizOutputEvent::dump(const std::string ind) const {
  cout << ind << "Viz output event:\n";
  std::string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
  cout << ind2 << "viz_mode: \t\t" << viz_mode << " [viz_mode_t] \t\t\n";
  cout << ind2 << "file_prefix_name: \t\t" << file_prefix_name << " [const char*] \t\t\n";
}


void VizOutputEvent::step() {
  switch (viz_mode) {
  case NO_VIZ_MODE:
    break;
  case ASCII_MODE:
    output_ascii_molecules();
    break;
  case CELLBLENDER_MODE:
    output_cellblender_molecules();
    break;
  default:
    assert(false);
  }
}


string VizOutputEvent::iterations_to_string(const uint64_t current_iterations, const uint64_t total_iterations) {
  stringstream res;

  uint64_t lli = 10;
  int ndigits;
  for (ndigits = 1; lli <= total_iterations && ndigits < 20; ndigits++) {
    lli *= 10;
  }

  res << setfill('0') << std::setw(ndigits) << current_iterations;
  return res.str();
}


FILE* VizOutputEvent::create_and_open_output_file_name() {

  long long current_iteration = round(event_time); // NOTE: usage of round might be a little shaky here, maybe we will need a better way how to get the iteration index
  //fprintf(stderr, "***dumps: %lld\n", current_iteration);
  const char* type_name = (viz_mode == ASCII_MODE) ? "ascii" : "cellbin";
  char* cf_name = CHECKED_SPRINTF(
      "%s.%s.%s.dat",
      file_prefix_name.c_str(),
      type_name,
      iterations_to_string(world->stats.get_current_iteration(), world->total_iterations).c_str()
   );
  assert(cf_name != nullptr);

  if (::make_parent_dir(cf_name)) {
    mcell_error(
        "Failed to create parent directory for CELLBLENDER-mode VIZ output.");
  }
  FILE *custom_file = ::open_file(cf_name, (viz_mode == ASCII_MODE) ? "w" : "wb");
  if (custom_file == nullptr)
    mcell_die();
  else {
    no_printf("Writing to file %s\n", cf_name);
  }
  free(cf_name);
  return custom_file;
}


void VizOutputEvent::compute_where_and_norm(
    const Partition& p, const Molecule& m,
    Vec3& where, Vec3& norm
) {
  const BNG::Species& species = world->get_all_species().get(m.species_id);

  // TODO: replace with new flags
  if (species.is_vol()) {
    // neither surface nor on grid
    where = m.v.pos;
    norm = Vec3(0);
  }
  else if (species.is_surf()) {
    const Wall& wall = p.get_wall(m.s.wall_index);
    const Vec3& wall_vert0 = p.get_geometry_vertex(wall.vertex_indices[0]);
    where = GeometryUtil::uv2xyz(m.s.pos, wall, wall_vert0);

    norm = Vec3(m.s.orientation) * wall.normal;
  }
  else {
    assert(false && "Unexpected molecule type");
    where = Vec3(0); // to silence compiler warnings
    norm = Vec3(0);
  }

  where *= world->config.length_unit;
}


void VizOutputEvent::output_ascii_molecules() {
  // assuming that fdlp->type == ALL_MOL_DATA
  FILE *custom_file = create_and_open_output_file_name();

  // simply go through all partitions and dump all molecules
  for (Partition& p: world->get_partitions()) {
    for (const Molecule& m: p.get_molecules()) {
      if (m.is_defunct()) {
        continue;
      }
      if (species_ids_to_visualize.count(m.species_id) == 0) {
        continue;
      }


      Vec3 where;
      Vec3 norm;
      compute_where_and_norm(p, m, where, norm);

      const BNG::Species& species = world->get_all_species().get(m.species_id);

      glm::dvec3 dwhere = where;
      glm::dvec3 dnorm = norm;
      assert(sizeof(dwhere.x) == sizeof(double));

      errno = 0;
      fprintf(custom_file, "%s %u %.9g %.9g %.9g %.9g %.9g %.9g\n",
          species.name.c_str(), m.id,
          dwhere.x, dwhere.y, dwhere.z,
          dnorm.x, dnorm.y, dnorm.z
      );
      assert(errno == 0);
    }
  }

  errno = 0;
  fclose(custom_file);
  assert(errno == 0);
}


void VizOutputEvent::output_cellblender_molecules() {
  // assuming that fdlp->type == ALL_MOL_DATA
  FILE *custom_file = create_and_open_output_file_name();
  float_t length_unit = world->config.length_unit;

  // sort all molecules by species
  uint species_count = world->get_all_species().get_count();

  // FIXME: rather change to partition and molecule indices
  typedef pair<const Partition*, const Molecule*> partition_molecule_ptr_pair_t;
  vector< vector<partition_molecule_ptr_pair_t> > volume_molecules_by_species;
  volume_molecules_by_species.resize(species_count);

  for (Partition& p: world->get_partitions()) {
    for (const Molecule& m: p.get_molecules()) {
      if (m.is_defunct()) {
        continue;
      }
      if (species_ids_to_visualize.count(m.species_id) == 0) {
        continue;
      }
      volume_molecules_by_species[m.species_id].push_back(partition_molecule_ptr_pair_t(&p, &m));
    }
  }

  /* Write file header */
  assert(sizeof(u_int) == sizeof(uint));
  uint cellbin_version = 1;
  fwrite(&cellbin_version, sizeof(cellbin_version), 1, custom_file);

  /* Write all the molecules whether EXTERNAL_SPECIES or not (for now) */
  for (species_id_t species_idx = 0; species_idx < world->get_all_species().get_count(); species_idx++) {
    // count of molecules for this species
    vector<partition_molecule_ptr_pair_t>& species_molecules = volume_molecules_by_species[species_idx];
    if (species_molecules.empty()) {
      continue;
    }

    /* Write species name: */
    const BNG::Species& species = world->get_all_species().get(species_idx);
    string mol_name = species.name;
    byte name_len = mol_name.length();
     fwrite(&name_len, sizeof(byte), 1, custom_file);
    fwrite(mol_name.c_str(), sizeof(char), name_len, custom_file);

     /* Write species type: */
    byte species_type = 0;
    if (species.is_surf()) {
      species_type = 1;
    }
    fwrite(&species_type, sizeof(species_type), 1, custom_file);

    /* write number of x,y,z floats for mol positions to follow: */
    uint n_floats = 3 * species_molecules.size();
    fwrite(&n_floats, sizeof(n_floats), 1, custom_file);

     /* Write positions of volume and surface molecules: */
     for (const partition_molecule_ptr_pair_t& partition_molecule_ptr_pair :species_molecules) {

       assert(partition_molecule_ptr_pair.second->is_vol() && "TODO - dump norm for surface molecules");

       Vec3 where;
       Vec3 norm;
       compute_where_and_norm(
           *partition_molecule_ptr_pair.first, *partition_molecule_ptr_pair.second,
           where, norm
       );

       // the values are always stored as double
       glm::dvec3 dwhere = where;
       assert(sizeof(dwhere.x) == sizeof(double));

       fwrite(&dwhere.x, sizeof(double), 1, custom_file);
       fwrite(&dwhere.y, sizeof(double), 1, custom_file);
       fwrite(&dwhere.z, sizeof(double), 1, custom_file);

       // TODO: store norm
     }

   }

  fclose(custom_file);
}


bool VizOutputEvent::should_visualize_all_species() const {

  // we are counting only specific vol & surf species, not reactive surfaces
  uint vol_surf_species_count = 0;
  for (const BNG::Species& s: world->get_all_species().get_species_vector()) {
    if (!is_species_superclass(s.name) && !s.is_reactive_surface()) {
      vol_surf_species_count++;
    }
  }

  return species_ids_to_visualize.size() == vol_surf_species_count;
}


void VizOutputEvent::to_data_model(Json::Value& mcell_node) const {
  // there can be just one viz_output section
  CONVERSION_CHECK(!mcell_node.isMember(KEY_VIZ_OUTPUT), "Only one viz_output section is allowed");

  Json::Value& viz_output = mcell_node[KEY_VIZ_OUTPUT];
  DMUtil::add_version(viz_output, VER_DM_2014_10_24_1638);

  viz_output[KEY_EXPORT_ALL] = should_visualize_all_species();
  viz_output[KEY_START] = DMUtil::f_to_string(event_time);
  viz_output[KEY_ALL_ITERATIONS] = cmp_eq(periodicity_interval, 1.0);
  viz_output[KEY_STEP] = DMUtil::f_to_string(periodicity_interval);
  viz_output[KEY_END] = to_string(world->total_iterations);

}

} // namespace mcell
