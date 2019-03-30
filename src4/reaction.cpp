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
#include "reaction.h"

using namespace std;

namespace mcell {

void reaction_t::dump(const string ind) const {
	cout << ind << "name: \t\t" << name << " [string] \t\t\n";
	cout << ind << "rate_constant: \t\t" << rate_constant << " [float_t] \t\t\n";
	cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";
}

} // namespace mcell
