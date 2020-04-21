/******************************************************************************
 *
 * Copyright (C) 2020 by
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

/**
 * This header includes all C++ classes of MCell API.
 * File should be used only from outside of this library to avoid cyclic includes.
 */

#ifndef API_MCELL_H
#define API_MCELL_H

#include "../api/common.h"


// data classes
#include "../api/component_type.h"
#include "../api/component_instance.h"
#include "../api/molecule_type.h"
#include "../api/molecule_instance.h"

#include "../api/geometry_object.h"
#include "../api/release_site.h"
#include "../api/species.h"


// classes with methods
#include "../api/model.h"
#include "../api/subsystem.h"
#include "../api/instantiation_data.h"

#endif // API_MCELL_H
