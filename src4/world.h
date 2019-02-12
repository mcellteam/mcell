/*
 * world.h
 *
 *  Created on: Feb 10, 2019
 *      Author: adam
 */

#ifndef SRC4_WORLD_H_
#define SRC4_WORLD_H_

#include <vector>

#include "partition.h"
#include "scheduler.h"
#include "species.h"

namespace mcell {

class world_t {
public:
  std::vector<partition_t> partitions;

  scheduler_t scheduler;

  std::vector<species_t> species;
};

} /* namespace mcell4 */

#endif /* SRC4_WORLD_H_ */
