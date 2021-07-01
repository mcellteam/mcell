/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

// this is an auxiliary file that provides std::vector functionality, used in viz_output.c

#include <vector>
#include <algorithm>

#include "c_vector.h"
#include "mcell_structs.h"


using namespace std;


typedef vector<void*> stl_vector_t;

c_vector_t* vector_create() {
  return new stl_vector_t();
}

void vector_push_back(c_vector_t* v, void* a) {
  ((stl_vector_t*)v)->push_back(a);
}

unsigned long long vector_get_size(c_vector_t* v) {
  return ((stl_vector_t*)v)->size();
}

void* vector_at(c_vector_t* v, unsigned long long i) {
  return ((stl_vector_t*)v)->at(i);
}

void vector_delete(c_vector_t* v) {
  delete (stl_vector_t*)v;
}

void vector_sort_by_mol_id(c_vector_t* v) {
  stl_vector_t& vec = *(stl_vector_t*)v;

  sort(vec.begin(), vec.end(),
      [](const void* a, const void* b) -> bool
  {
      return ((abstract_molecule*)a)->id < ((abstract_molecule*)b)->id;
  });
}

