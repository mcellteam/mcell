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
#include "release_event.h"

extern "C" {
#include "logging.h"
}

// checking major converstion blocks
#define CHECK(a) do { if(!(a)) return false; } while (0)

// checking assumptions
#define CHECK_PROPERTY(cond) do { if(!(cond)) { mcell_log_conv_error("Expected '%s' is false. (%s:%d)", #cond, __FILE__, __LINE__); return false; } } while (0)

// asserts - things that can never occur and will 'never' be supported


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


static mat4x4 t_matrix_to_mat4x4(const double src[4][4]) {
	mat4x4 res;

	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			res[x][y] = src[x][y];
		}
	}

	return res;
}



bool mcell3_world_converter::convert(volume* s) {

	world = new world_t();

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

		mcell3_species_id_map[new_species.mcell3_species_id] = new_species.species_id;
  }

	return true;
}

bool mcell3_world_converter::convert_release_events(volume* s) {

	// -- schedule_helper -- (as volume.releaser)
	schedule_helper* releaser = s->releaser;

	CHECK_PROPERTY(releaser->next_scale == nullptr);
	CHECK_PROPERTY(releaser->dt == 1);
	CHECK_PROPERTY(releaser->dt_1 == 1);
	CHECK_PROPERTY(releaser->now == 0);
	CHECK_PROPERTY(releaser->count == 1);
	CHECK_PROPERTY(releaser->index == 0);

  for (int i = -1; i < releaser->buf_len; i++) {
  	int k = 0;
    for (abstract_element *aep = (i < 0) ? releaser->current
                                                : releaser->circ_buf_head[i];
         aep != NULL; aep = aep->next) {

    	release_event_t* event = new release_event_t();

    	// -- release_event_queue --
			release_event_queue *req = (release_event_queue *)aep;


			event->event_time = req->event_time;

			// -- release_site --
			release_site_obj* rel_site = req->release_site;

			assert(rel_site->location != nullptr);
			event->location = *rel_site->location;
			event->species_id = get_new_species_id(rel_site->mol_type->species_id);

			CHECK_PROPERTY(rel_site->release_number_method == 0);
			CHECK_PROPERTY(rel_site->release_shape == 0);
			CHECK_PROPERTY(rel_site->orientation == 0);

			event->release_number = rel_site->release_number;

			CHECK_PROPERTY(rel_site->mean_diameter == 0); // temporary
			CHECK_PROPERTY(rel_site->concentration == 0); // temporary
			CHECK_PROPERTY(rel_site->standard_deviation == 0); // temporary
			assert(rel_site->diameter != nullptr);
			CHECK_PROPERTY(vec3_t(*rel_site->diameter) == vec3_t(0,0,0)); // temporary
			CHECK_PROPERTY(rel_site->region_data == nullptr); // temporary?
			CHECK_PROPERTY(rel_site->mol_list == nullptr);
			CHECK_PROPERTY(rel_site->release_prob == 1); // temporary
			// rel_site->periodic_box - ignoring?
			// rel_site->pattern - TODO - is not null
			event->name = rel_site->name;
			// rel_site->graph_pattern - TODO - is not null - NFSim?

    	// -- release_event_queue -- (again)
			CHECK_PROPERTY(t_matrix_to_mat4x4(req->t_matrix) == mat4x4(1) && "only identity matrix for now");
			CHECK_PROPERTY(req->train_counter == 0);
	    CHECK_PROPERTY(req->train_high_time == 0);
			CHECK_PROPERTY(req->next == nullptr);

			world->scheduler.schedule_event(event);
    }
  }

  // -- schedule_helper -- (again)
  CHECK_PROPERTY(releaser->current_count ==	0);
  CHECK_PROPERTY(releaser->current == nullptr);
  CHECK_PROPERTY(releaser->current_tail == nullptr);
  CHECK_PROPERTY(releaser->defunct_count == 0);
  CHECK_PROPERTY(releaser->error == 0);
  CHECK_PROPERTY(releaser->depth == 0);

	return true;
}


} /* namespace mcell */
