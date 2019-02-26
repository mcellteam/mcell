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

extern "C" {
#include "logging.h"
#include "mem_util.h"
#include "util.h" // MCell 3
#include "mcell_structs.h"
}

#include "viz_output_event.h"
#include "world.h"


using namespace std;

namespace mcell {


void viz_output_event_t::dump(const std::string indent) {
	cout << indent << "Viz output event:\n";
	std::string ind2 = indent + "  ";
	base_event_t::dump(ind2);
	cout << ind2 << "viz_mode: \t\t" << viz_mode << " [viz_mode_t] \t\t\n";
	cout << ind2 << "file_prefix_name: \t\t" << file_prefix_name << " [const char*] \t\t\n";
}


void viz_output_event_t::step() {
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

FILE* viz_output_event_t::create_and_open_output_file_name() {
	int ndigits = digits_for_file_suffix(world->iterations);
	long long current_iteration = round(event_time / world->time_unit); // FIXME: usage of round might be a little shaky here, maybe we will need a better way how to get the iteration index
  //fprintf(stderr, "***dumps: %lld\n", current_iteration);
	const char* type_name = (viz_mode == ASCII_MODE) ? "ascii" : "cellbin";
  char* cf_name =
  		CHECKED_SPRINTF("%s.%s.%.*lld.dat", file_prefix_name, type_name, ndigits, current_iteration);
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


void viz_output_event_t::output_ascii_molecules() {
	// assuming that fdlp->type == ALL_MOL_DATA
	FILE *custom_file = create_and_open_output_file_name();

	// simply go through all partitions and dump all molecules
	uint64_t id = 0;
	for (partition_t& p: world->partitions) {
		for (volume_molecule_t& m: p.volume_molecules) {
			std::string species_name = world->species[m.species_id].name;
#if FLOAT_T_BYTES == 8
			// TODO: norm
			errno = 0;
      fprintf(custom_file, "%s %lu %.9g %.9g %.9g %.9g %.9g %.9g\n",
      		species_name.c_str(), id, m.pos.x, m.pos.y, m.pos.z, 0.0, 0.0, 0.0
			);
      assert(errno == 0);
#else
#error "Marker for float type"
#endif
			id++;
		}
	}

	errno = 0;
  fclose(custom_file);
  assert(errno == 0);
}

void viz_output_event_t::output_cellblender_molecules() {
	// assuming that fdlp->type == ALL_MOL_DATA
	FILE *custom_file = create_and_open_output_file_name();

	// sort all molecules by species
	uint32_t species_count = world->species.size();

	vector< vector<volume_molecule_t*> > volume_molecules_by_species;
	volume_molecules_by_species.resize(species_count);

	for (partition_t& p: world->partitions) {
		for (volume_molecule_t& m: p.volume_molecules) {
			volume_molecules_by_species[m.species_id].push_back(&m);
		}
	}

  /* Write file header */
	assert(sizeof(u_int) == sizeof(uint32_t));
	uint32_t cellbin_version = 1;
  fwrite(&cellbin_version, sizeof(cellbin_version), 1, custom_file);

  /* Write all the molecules whether EXTERNAL_SPECIES or not (for now) */
  for (species_id_t species_idx = 0; species_idx < world->species.size(); species_idx++) {
  	// count of molecules for this species
  	vector<volume_molecule_t*>& species_molecules = volume_molecules_by_species[species_idx];
  	if (species_molecules.empty()) {
  		continue;
  	}

    /* Write species name: */
    string mol_name = world->species[species_idx].name;
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
		uint32_t n_floats = 3 * species_molecules.size();
    fwrite(&n_floats, sizeof(n_floats), 1, custom_file);

     /* Write positions of volume and surface molecules: */
     float pos_x = 0.0;
     float pos_y = 0.0;
     float pos_z = 0.0;
     for (volume_molecule_t* mp : species_molecules) {
    	 // TODO: many specific variants missing
			 pos_x = mp->pos.x;
			 pos_y = mp->pos.y;
			 pos_z = mp->pos.z;

			 // this scould be already incorporated
       /*pos_x *= world->length_unit;
       pos_y *= world->length_unit;
       pos_z *= world->length_unit;*/

       fwrite(&pos_x, sizeof(pos_x), 1, custom_file);
       fwrite(&pos_y, sizeof(pos_y), 1, custom_file);
       fwrite(&pos_z, sizeof(pos_z), 1, custom_file);
     }
   }

  fclose(custom_file);

}

} /* namespace mcell */
