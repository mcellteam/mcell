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

%module pymcell3_4 

%import stdint.i

typedef unsigned int uint;

%{

#include "world.h"
#include "mcell3_world_converter.h" 

using namespace MCell;

/* This function matches the prototype of the normal C callback
   function for our widget. However, we use the clientdata pointer
   for holding a reference to a Python callable object. */

// TODO: where should I put these callbacks?

#include "callback_info.h"
%}

// TODO: why does simple import not work?
//%import "defines.h"
typedef uint partition_index_t;
typedef uint vertex_index_t;

%import "callback_structs.h"


// Grab a Python function object as a Python object.
%typemap(in) PyObject *pyfunc {
  if (!PyCallable_Check($source)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $target = $source;
}

// Type mapping for grabbing a FILE * from Python
%typemap(in) FILE * {
  if (!PyFile_Check($source)) {
      PyErr_SetString(PyExc_TypeError, "Need a file!");
      return NULL;
  }
  $target = PyFile_AsFile($source);
}

namespace MCell {

typedef double float_t;

struct vec3_t {
  vec3_t();
  vec3_t(const float_t x_, const float_t y_, const float_t z_);

  float_t x, y, z; # should be float_t
};

class Partition {
public:
  Partition(
      const vec3_t origin_,
      const WorldConstants& world_constants_,
      SimulationStats& simulation_stats_
  );
  
  int get_geometry_vertex_count();
  vec3_t& get_geometry_vertex(vertex_index_t i);
  
  void add_vertex_move(vertex_index_t vertex_index, const vec3_t& translation_vec);
  void apply_vertex_moves();
};

class WorldConstants {
public:
  float_t length_unit;
 
};
  
class World {
public:
  void run_simulation(const bool dump_initial_state = false);
  void run_n_iterations(const uint64_t num_iterations, const uint64_t output_frequency);
  void end_simulation();
  
  void register_wall_hit_callback_internal(wall_hit_callback_func func, void* clientdata_);
  
  WorldConstants& get_world_constants();
  
  Partition& get_partition(partition_index_t i);
  
  void dump();
};

class MCell3WorldConverter {
public:
  bool convert(MCELL_STATE* s);
  World* world;
};

} // namespace MCell

// Attach a new method to our plot widget for adding Python functions
%extend World {
   // Set a Python function object as a callback function
   // Note : PyObject *pyfunc is remapped with a typempap
   void register_wall_hit_callback(PyObject *pyfunc) {
     self->register_wall_hit_callback_internal(py_callback_wall_hit, (void *) pyfunc);
     Py_INCREF(pyfunc);
   }
}
