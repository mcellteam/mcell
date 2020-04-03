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

#ifndef SRC4_MOL_COUNT_EVENT_H_
#define SRC4_MOL_COUNT_EVENT_H_

#include "base_event.h"

#include "count_buffer.h"

namespace MCell {

class Partition;
class Molecule;

enum class CountType {
  Invalid,
  World,
  EnclosedInObject
};

class MolCountInfo {
public:
  MolCountInfo(const count_buffer_id_t buffer_id_)
    : type(CountType::Invalid),
      buffer_id(buffer_id_),
      orientation(ORIENTATION_NOT_SET),
      species_id(SPECIES_ID_INVALID),
      geometry_object_id(GEOMETRY_OBJECT_ID_INVALID)
      {
    assert(buffer_id != COUNT_BUFFER_ID_INVALID);
  }

  void dump(const std::string ind);

  CountType type;

  // count buffer objects are owned by World
  count_buffer_id_t buffer_id;

  orientation_t orientation;
  species_id_t species_id;

  // valid when type is EnclosedInObject
  geometry_object_id_t geometry_object_id;
};


/**
 * Dumps counts of molecules.
 */
class MolCountEvent: public BaseEvent {
public:
  MolCountEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_MOL_COUNT),
      world(world_) {
  }
  virtual ~MolCountEvent() {}

  virtual void step();

  // FIMXE: all events' dumps must be const and should use just 'ind' as arg name
  virtual void dump(const std::string ind);

  void add_mol_count_info(const MolCountInfo& info) {
    mol_count_infos.push_back(info);
  }

  std::vector<MolCountInfo> mol_count_infos;

  World* world;
};

} // namespace mcell

#endif // SRC4_MOL_COUNT_EVENT_H_
