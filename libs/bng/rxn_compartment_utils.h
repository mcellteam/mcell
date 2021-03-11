/*
 * rxn_compartment_utils.h
 *
 *  Created on: Oct 27, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_RXN_COMPARTMENT_UTILS_H_
#define LIBS_BNG_RXN_COMPARTMENT_UTILS_H_


namespace BNG {
class BNGData;
class RxnRule;

// - MCell-specific, replaces compartments (including @IN and @OUT) with appropriate orientations
// - keeps the original compartment information there but it is ignored
//   (considered to be COMPARTMENT_ID_NONE)
// - the conversion is equivalent as long as there are no reactions that change the MCell
//   orientation of surface molecules, the assumption is that all surface molecules are oriented
//   outwards (up) all the time
// - returns nonempty string if there was and error
std::string process_compartments_and_set_orientations(const BNGData& bng_data, RxnRule& r);

}

#endif /* LIBS_BNG_RXN_COMPARTMENT_UTILS_H_ */
