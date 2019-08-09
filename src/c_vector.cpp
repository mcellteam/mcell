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

