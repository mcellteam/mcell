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
#include <stdio.h>
#include <errno.h>

//extern "C" {
#include "logging.h"
#include "mem_util.h"
#include "util.h"
#include "mcell_structs.h"
//}

#include "viz_output_event.h"
#include "geometry.h"
#include "world.h"


#include "geometry_utils.inc"

using namespace std;

namespace MCell {

void VizOutputEvent::dump(const std::string indent) {
  cout << indent << "Viz output event:\n";
  std::string ind2 = indent + "  ";
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


static int digits_for_file_suffix(uint64_t iterations) {
  uint64_t lli = 10;
  int ndigits;
  for (ndigits = 1; lli <= iterations && ndigits < 20; ndigits++) {
    lli *= 10;
  }
  return ndigits;
}


FILE* VizOutputEvent::create_and_open_output_file_name() {
  int ndigits = digits_for_file_suffix(world->iterations);
  long long current_iteration = round(event_time); // NOTE: usage of round might be a little shaky here, maybe we will need a better way how to get the iteration index
  //fprintf(stderr, "***dumps: %lld\n", current_iteration);
  const char* type_name = (viz_mode == ASCII_MODE) ? "ascii" : "cellbin";
  char* cf_name =
      CHECKED_SPRINTF("4%s.%s.%.*lld.dat", file_prefix_name, type_name, ndigits, current_iteration);
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
    vec3_t& where, vec3_t& norm
) {
  const Species& species = world->get_species(m.species_id);

  if ((species.flags & NOT_FREE) == 0) {
    // neither surface nor on grid
    where = m.v.pos;
    norm = vec3_t(0);
  }
  else if ((species.flags & ON_GRID) != 0) {
    const Wall& wall = p.get_wall(m.s.wall_index);
    const vec3_t& wall_vert0 = p.get_geometry_vertex(wall.vertex_indices[0]);
    where = GeometryUtil::uv2xyz(m.s.pos, wall, wall_vert0);

    norm = vec3_t(m.s.orientation) * wall.normal;
  }
  else {
    assert(false && "Unexpected molecule type");
    where = vec3_t(0); // to silence compiler warnings
    norm = vec3_t(0);
  }

  where *= world->world_constants.length_unit;
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


      vec3_t where;
      vec3_t norm;
      compute_where_and_norm(p, m, where, norm);

      const Species& species = world->get_species(m.species_id);

#if FLOAT_T_BYTES == 8
      // TODO: norm
      // also some test that the output is correct is needed
      errno = 0;
      fprintf(custom_file, "%s %u %.9g %.9g %.9g %.9g %.9g %.9g\n",
          species.name.c_str(), m.id,
          where.x, where.y, where.z,
          norm.x, norm.y, norm.z
      );
      assert(errno == 0);
#else
#error "Marker for float type"
#endif
    }
  }

  errno = 0;
  fclose(custom_file);
  assert(errno == 0);
}


void VizOutputEvent::output_cellblender_molecules() {
  // assuming that fdlp->type == ALL_MOL_DATA
  FILE *custom_file = create_and_open_output_file_name();
  float_t length_unit = world->world_constants.length_unit;

  // sort all molecules by species
  uint species_count = world->get_species().size();

  // FIXME: rather change to partition and molecule indices
  typedef pair<const Partition*, const Molecule*> partition_molecule_ptr_pair_t;
  vector< vector<partition_molecule_ptr_pair_t> > volume_molecules_by_species;
  volume_molecules_by_species.resize(species_count);

  for (Partition& p: world->get_partitions()) {
    for (const Molecule& m: p.get_molecules()) {
      if (m.is_defunct()) {
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
  for (species_id_t species_idx = 0; species_idx < world->get_species().size(); species_idx++) {
    // count of molecules for this species
    vector<partition_molecule_ptr_pair_t>& species_molecules = volume_molecules_by_species[species_idx];
    if (species_molecules.empty()) {
      continue;
    }

    /* Write species name: */
    string mol_name = world->get_species(species_idx).name;
    byte name_len = mol_name.length();
     fwrite(&name_len, sizeof(byte), 1, custom_file);
    fwrite(mol_name.c_str(), sizeof(char), name_len, custom_file);

     /* Write species type: */
    byte species_type = 0;
    /*TODO: if ((amp->properties->flags & ON_GRID) != 0) {
      species_type = 1;
    }*/
    fwrite(&species_type, sizeof(species_type), 1, custom_file);

    /* write number of x,y,z floats for mol positions to follow: */
    uint n_floats = 3 * species_molecules.size();
    fwrite(&n_floats, sizeof(n_floats), 1, custom_file);

     /* Write positions of volume and surface molecules: */
     for (const partition_molecule_ptr_pair_t& partition_molecule_ptr_pair :species_molecules) {

       assert(partition_molecule_ptr_pair.second->is_vol() && "TODO - dump norm for surface molecules");

       vec3_t where;
       vec3_t norm;
       compute_where_and_norm(
           *partition_molecule_ptr_pair.first, *partition_molecule_ptr_pair.second,
           where, norm
       );

#if FLOAT_T_BYTES == 8
       fwrite(&where.x, sizeof(float_t), 1, custom_file);
       fwrite(&where.y, sizeof(float_t), 1, custom_file);
       fwrite(&where.z, sizeof(float_t), 1, custom_file);
#else
#error "Marker for float type"
#endif
     }

   }

  fclose(custom_file);
}

} // namespace mcell