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

// this is an auxiliary file that provides std::vector functionality that contains pointers as elements
// used in viz_output.c

#ifndef C_VECTOR_H
#define C_VECTOR_H

#ifdef __cplusplus
//extern "C" {
#endif

typedef void c_vector_t;

c_vector_t* vector_create();
void vector_push_back(c_vector_t* v, void* a);
unsigned long long vector_get_size(c_vector_t* v);
void* vector_at(c_vector_t* v, unsigned long long i);
void vector_delete(c_vector_t* v);


void vector_sort_by_mol_id(c_vector_t* v);


#ifdef __cplusplus
//}
#endif

#endif // C_VECTOR_H

