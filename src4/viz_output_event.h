/*
 * viz_output_event.h
 *
 *  Created on: Feb 19, 2019
 *      Author: adam
 */

#ifndef SRC4_VIZ_OUTPUT_EVENT_H_
#define SRC4_VIZ_OUTPUT_EVENT_H_

#include "base_event.h"

namespace mcell {

class viz_output_event: public base_event_t {
public:
	viz_output_event()
		: base_event_t(EVENT_TYPE_VIZ_OUTPUT) {
	}
	virtual void step();
	virtual void dump() { /*empty*/ }
};

} /* namespace mcell */

#endif /* SRC4_VIZ_OUTPUT_EVENT_H_ */
