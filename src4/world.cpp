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

#include <time.h>
#include <sys/time.h> // Linux include
#include <sys/resource.h> // Linux include
#include <fenv.h> // Linux include

extern "C" {
#include "rng.h" // MCell 3
#include "logging.h"
}

#include "world.h"
#include "end_simulation_event.h"
#include "defragmentation_event.h"

using namespace std;

const double USEC_IN_SEC = 1000000.0;

namespace mcell {

world_t::world_t()
  : current_iteration(0),
    iterations(0),
    seed_seq(0),
    next_wall_id(0),
    next_geometry_object_id(0)
{
  // TODO: initialize rest of members
  world_constants.partition_edge_length = PARTITION_EDGE_LENGTH_DEFAULT;
  world_constants.subpartitions_per_partition_dimension = SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT;
}


void world_t::init_fpu() {
#ifdef NDEBUG
  // we do not want to be making extra checks for division by zero
  // all places where such a case can occur is marked with comment POSSIBLE ZERO DIV
  fedisableexcept(FE_DIVBYZERO);

  static float_t a = 1;
  static float_t b = 0;
  if (a/b != INFINITY) {
    mcell_error("Error: division by zero is expected to return INFINITY but does not!");
    exit(1);
  }
#endif
}

void world_t::init_world_constants() {

  if (species.empty()) {
    // for developers: init_world_constants must be called after species and reactions conversion
    mcell_log("Error: there must be at lease one species!");
    exit(1);
  }

  // create map for fast reaction searches
  for (reaction_t& r: reactions) {
    assert(r.reactants.size() == 1 || r.reactants.size() == 2); // only bimolecular reactions are supported now

    if (r.reactants.size() == 1) {
      // for now we only support only one outcome of a bimolecular reaction
      assert(unimolecular_reactions_map.count(r.reactants[0].species_id) == 0);
      unimolecular_reactions_map[r.reactants[0].species_id] = &r;
    }
    else {
      // check - for now we only support only one outcome of a bimolecular reaction
      if (bimolecular_reactions_map.count(r.reactants[0].species_id) != 0) {
        assert(bimolecular_reactions_map[r.reactants[0].species_id].count(r.reactants[1].species_id) == 0);
      }
      bimolecular_reactions_map[r.reactants[0].species_id][r.reactants[1].species_id] = &r;
      bimolecular_reactions_map[r.reactants[1].species_id][r.reactants[0].species_id] = &r;
    }
  }

  // just to make sure that we have an item for all the species
  for (species_t& s: species) {
    bimolecular_reactions_map.insert( std::make_pair(s.species_id, species_reaction_map_t()) );
  }
  assert(bimolecular_reactions_map.size() == species.size());

  world_constants.init(&unimolecular_reactions_map, &bimolecular_reactions_map);
}


void world_t::init_simulation() {

  init_fpu();

  cout << "Partitions contain " <<  world_constants.subpartitions_per_partition_dimension << "^3 subvolumes.";
  assert(partitions.size() == 1 && "Initial parition must have been created, only 1 for now");
}


static uint64_t determine_output_frequency(uint64_t iterations) {
  uint64_t frequency;

  if (iterations < 10)
    frequency = 1;
  else if (iterations < 1000)
    frequency = 10;
  else if (iterations < 100000)
    frequency = 100;
  else if (iterations < 10000000)
    frequency = 1000;
  else if (iterations < 1000000000)
    frequency = 10000;
  else
    frequency = 100000;

  return frequency;
}


static double tousecs(timeval& t) {
  return (double)t.tv_sec * USEC_IN_SEC + (double)t.tv_usec;
}


static double tosecs(timeval& t) {
  return (double)t.tv_sec + (double)t.tv_usec/USEC_IN_SEC;
}


bool world_t::run_simulation() {

  init_simulation(); // must be the first one

  dump();

  // create defragmentation events
  defragmentation_event_t* defragmentation_event = new defragmentation_event_t(this);
  defragmentation_event->event_time = DEFRAGMENTATION_PERIODICITY;
  defragmentation_event->periodicity_interval = DEFRAGMENTATION_PERIODICITY;
  scheduler.schedule_event(defragmentation_event);

  // create event that will terminate our simulation
  end_simulation_event_t* end_event = new end_simulation_event_t();
  end_event->event_time = iterations;
  scheduler.schedule_event(end_event);

  bool end_simulation = false;
  float_t time = TIME_SIMULATION_START;
  float_t previous_time;
  current_iteration = 0;
  uint32_t output_frequency = determine_output_frequency(iterations);
  timeval last_timing_time = {0, 0};

  cout << "Iterations: " << current_iteration << " of " << iterations << "\n";

  rusage sim_start_time;
  getrusage(RUSAGE_SELF, &sim_start_time);

  do {
    previous_time = time;

    // this is where events get executed
    time = scheduler.handle_next_event(end_simulation);

    // report progress
    if (time > previous_time) {
      current_iteration++;
      if (current_iteration % output_frequency == 0) {

        cout << "Iterations: " << current_iteration << " of " << iterations;

        timeval curr_timing_time;
        gettimeofday(&curr_timing_time, NULL);
        if (last_timing_time.tv_usec > 0) {
          double time_diff = tousecs(curr_timing_time) - tousecs(last_timing_time);
          time_diff /= (double)output_frequency;
          cout << " (" << 1000000.0/time_diff << " iter/sec)";
        }
        last_timing_time = curr_timing_time;

        cout << "\n";
      }
    }

  } while (!end_simulation);

  cout << "Iteration " << current_iteration << ", simulation finished successfully\n";

  // report final time
  rusage run_time;
  getrusage(RUSAGE_SELF, &run_time);
  cout << "Simulation CPU time = "
    << tosecs(run_time.ru_utime) - tosecs(sim_start_time.ru_utime) <<  "(user) and "
    << tosecs(run_time.ru_stime) - tosecs(sim_start_time.ru_stime) <<  "(system)\n";

  return true;
}


partition_index_t world_t::get_partition_index_for_pos(const vec3_t& pos) {
  // for now very simply search all partitions
  // we expect that parittion boundaries are precise
  partition_index_t found_index = PARTITION_INDEX_INVALID;
  for (partition_index_t i = 0; i < partitions.size(); i++) {
    const partition_t& p = partitions[i];
    if ( glm::all( glm::greaterThanEqual(pos, p.get_origin_corner() )) &&
        glm::all( glm::lessThan(pos, p.get_opposite_corner()) ) ) {
      assert(found_index == PARTITION_INDEX_INVALID && "Single point can be in just one partition");
      found_index = i;
#ifdef NDEBUG
      break; // we found our partition
#endif
    }
  }
  return found_index;
}


partition_vertex_index_pair_t world_t::add_geometry_vertex(const vec3_t& pos) {
  // to which partition it belongs
  partition_index_t partition_index = get_partition_index_for_pos(pos);
  if (partition_index == PARTITION_INDEX_INVALID) {
    mcell_log("Error: only a single partition is supported for now, vertex %s is out of bounds", pos.to_string().c_str());
  }

  vertex_index_t vertex_index = partitions[partition_index].add_geometry_vertex(pos);

  return partition_vertex_index_pair_t(partition_index, vertex_index);
}


// adds a new geometry object with its walls, sets unique ids for the walls and objects
// note: there is a lot of potentially uncenecssary copying of walls, can be optimized
void world_t::add_geometry_object(
    const geometry_object_t& obj,
    const vector<wall_t>& walls, // the vertices for walls are contained in walls_vertices
    const vector<vector<partition_vertex_index_pair_t>>& walls_vertices
) {
  assert(!walls.empty());
  assert(walls.size() == walls_vertices.size());

  partition_index_t partition_index = walls_vertices[0][0].first;
  partition_t& p = partitions[partition_index];

  geometry_object_t& new_obj = p.add_geometry_object(obj, next_geometry_object_id);

  for (uint32_t i = 0; i < walls.size(); i++) {
    wall_t new_wall = walls[i];

    // check that all vertices are in the same partition (the previous code assumes this)
    assert(walls_vertices[i].size() == VERTICES_IN_TRIANGLE);
    for (uint32_t k = 0; k < VERTICES_IN_TRIANGLE; k++) {
      vertex_index_t vertex_index = walls_vertices[i][k].second;

      // check that we fit int the parition
      if (walls_vertices[i][k].first != partition_index) {
        vec3_t pos = p.get_geometry_vertex(vertex_index);
        mcell_log("Error: only a single partition is supported for now, vertex %s is out of bounds", pos.to_string().c_str());
      }

      // set walls to the object
      new_wall.vertex_indices[k] = vertex_index;
    }

    // add wall to partition
    wall_index_t new_wall_index; // !!! FIXME: the ordering is wrong,  p.add_wall needs vertices!
    p.add_wall(new_wall, next_wall_id, new_wall_index);

    // set index of the contained wall
    new_obj.wall_indices.push_back(new_wall_index);
  }
}


void world_t::dump() {
  world_constants.dump();
  // species
  species_t::dump_array(species);

  // paritions
  for (partition_t& p: partitions) {
    p.dump();
  }
}

} // namespace mcell

