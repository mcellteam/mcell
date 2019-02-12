/*
 * event.h
 *
 *  Created on: Jan 30, 2019
 *      Author: adam
 */

#ifndef _EVENT_H_
#define _EVENT_H_

namespace mcell {

/**
 * Base class for all events.
 */
class base_event_t {
public:
	virtual ~base_event_t();
	virtual void step() = 0;
};

}

#endif /* V4SRC_EVENT_H_ */
