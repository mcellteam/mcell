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

#include "legacy_callback_info.h"
%}

// TODO: why does simple import not work?
//%import "defines.h"
typedef uint partition_id_t;
typedef uint vertex_index_t;
typedef uint molecule_id_t;
typedef uint geometry_object_id_t;
typedef uint wall_id_t;


struct rng_state {
  unsigned int randcnt;
  unsigned long long aa;
  unsigned long long bb;
  unsigned long long cc;
  unsigned long long randrsl[256];
  unsigned long long mm[256];
  unsigned long long rngblocks;
};

// FIXME: this does not work either
//%import "callback_structs.h"


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

namespace BNG {
  typedef double float_t;

  class SpeciesContainer { 
    Species* find_species_by_name(const char* name);
  };
      
  %extend SpeciesContainer {      
    SpeciesContainer() {
     // not really used
    } 
  }
      
  class BNGData {
  };
      
  class Species {  
  public:
    Species(const BNGData& data);
    void set_color(float_t r, float_t g, float_t b);
    void set_scale(float_t s);
  };
} // namespace BNG

namespace MCell {

typedef double float_t;

struct Vec3 {
  Vec3();
  Vec3(const float_t x_, const float_t y_, const float_t z_);

  float_t x, y, z; # should be float_t
};

namespace API {
struct WallHitInfo {
  molecule_id_t molecule_id;
  geometry_object_id_t geometry_object_id;
  wall_id_t wall_id;
  float_t time;
  Vec3 pos;
  float_t time_before_hit;
  Vec3 pos_before_hit;
};
}


class SimulationConfig {
public:
  float_t length_unit;
};

class SpeciesFlagsAnalyzer {
};

class Partition {
public:
  // ctor definition is needed, default variant is generated instead by swig 
  Partition(
      const partition_id_t id_,  
      const Vec3 origin_,
      const SimulationConfig& config_,
      BNG::BNGEngine& bng_engine_,
      SimulationStats& stats_,
      SpeciesFlagsAnalyzer& species_flags_analyzer_
  );
  
  int get_geometry_vertex_count();
  Vec3& get_geometry_vertex(vertex_index_t i);
  
  void add_vertex_move(vertex_index_t vertex_index, const Vec3& translation_vec);
  void apply_vertex_moves();
  
  void dump();
};


class World {
public:
  void run_simulation(const bool dump_initial_state = false);
  void run_n_iterations(const uint64_t num_iterations, const uint64_t output_frequency);
  void end_simulation();
  
  void register_wall_hit_callback_internal(wall_hit_callback_func func, void* clientdata_, const char* object_name);
  
  
  void enable_wall_hit_counting();
  uint get_wall_hit_array_size();
  const API::WallHitInfo& get_wall_hit_array_item(uint index);
  void clear_wall_hit_array();
  
  Partition& get_partition(partition_id_t i);
  
  void dump();
  
  BNG::SpeciesContainer& get_all_species();
  
  SimulationConfig config;
};

class MCell3WorldConverter {
public:
  bool convert(MCELL_STATE* s);
  World* world;
};


// Attach a new method to our plot widget for adding Python functions
%extend World {
   // Set a Python function object as a callback function
   // Note : PyObject *pyfunc is remapped with a typemap
   void register_wall_hit_callback(PyObject *pyfunc, const char* object_name = "") {
     self->register_wall_hit_callback_internal(py_callback_wall_hit, (void *) pyfunc, object_name);
     Py_INCREF(pyfunc);
   }

  void export_visualization_datamodel_to_dir(const char* prefix) {
  	self->export_data_model_to_dir(std::string(prefix), true);
  }

  void export_datamodel(const char* file) {
  	self->export_data_model(std::string(file), false);
  }
   
    // only in Python API
    BNG::Species* find_species_by_name(const char* name) {
     	return self->get_all_species().find_species_by_name(name);
    }
}

} // namespace MCell
