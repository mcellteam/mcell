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

#include "base_event.h"

using namespace std;

namespace mcell {

void base_event_t::dump(const std::string indent) {
	cout << indent << "event_time: \t\t" << event_time << " [float_t] \t\t\n";
	cout << indent << "periodicity_interval: \t\t" << periodicity_interval << " [float_t] \t\t\n";
	cout << indent << "type_index: \t\t" << type_index << " [event_type_index_t] \t\t\n";
}

} // namespace mcell
