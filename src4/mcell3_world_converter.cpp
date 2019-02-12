/*
 * mcell3_world_converter.cpp
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#include <stdarg.h>
#include <stdlib.h>

#include "mcell_structs.h"

#include "mcell3_world_converter.h"

#include "world.h"

extern "C" {
#include "logging.h"
}


#define CHECK(a) do { if(!(a)) return false; } while (0)

// holds global class
mcell::mcell3_world_converter g_converter;

bool convert_mcell3_volume_to_mcell3_world(volume* s) {
	return g_converter.convert(s);
}

void mcell_log_conv_warning(char const *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  std::string fmt_w_warning = std::string ("Conversion warning: ") + fmt;
  mcell_logv_raw(fmt_w_warning.c_str(), args);
  va_end(args);
}

void mcell_log_conv_error(char const *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  std::string fmt_w_warning = std::string ("Conversion error: ") + fmt;
  mcell_logv_raw(fmt_w_warning.c_str(), args);
  va_end(args);
}

namespace mcell {

static const char* get_sym_name(const sym_entry *s) {
	assert(s != nullptr);
	assert(s->name != nullptr);
	return s->name;
}



bool mcell3_world_converter::convert(volume* s) {
	// 1) species
	CHECK(convert_species(s));

  // 2) release events
	CHECK(convert_release_events(s));


	return true;
}


bool mcell3_world_converter::convert_species(volume* s) {

  for (int i = 0; i < s->n_species; i++) {
  	species* spec = s->species_list[i];
  	// not sure what to do with these superclasses
  	if (spec == s->all_mols || spec == s->all_surface_mols || spec == s->all_surface_mols)
  		continue;

  	if (spec->flags != 0) {
  		mcell_log_conv_error("Species with index %d: only flags == 0 are supported");
  		return false;
  	}

		species_t new_species;

		new_species.species_id = world->species.size(); // id corresponds to the index in the species array
		new_species.mcell3_species_id = spec->species_id;
		new_species.D = spec->D;
		new_species.name = get_sym_name(spec->sym);
		new_species.space_step = spec->space_step;
		new_species.time_step = spec->time_step;

		world->species.push_back(new_species);
  }

	return true;
}

bool mcell3_world_converter::convert_release_events(volume* s) {



	return true;
}


} /* namespace mcell */
