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

// this is an auxiliary file that provides std::vector functionality that contains pointers as elements
// used in viz_output.c

#ifndef C_VECTOR_H
#define C_VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

typedef void c_vector_t;

c_vector_t* vector_create();
void vector_push_back(c_vector_t* v, void* a);
unsigned long long vector_get_size(c_vector_t* v);
void* vector_at(c_vector_t* v, unsigned long long i);
void vector_delete(c_vector_t* v);


void vector_sort_by_mol_id(c_vector_t* v);


#ifdef __cplusplus
}
#endif

#endif // C_VECTOR_H

