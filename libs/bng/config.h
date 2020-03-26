/*
 * config.h
 *
 *  Created on: Mar 26, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_CONFIG_H_
#define LIBS_BNG_CONFIG_H_

// ---------------------------------- configurability ------------------------------
/**
 * The goal of this library of for it to be general enough so that different
 * simulation engines can use it.
 *
 * One such example is what information does the specific engine needs about
 * species. E.g. MCell species uses variables space_step and time_step.
 * These are most probably of little use.
 *
 * One approach to configurability are templates. However the issue with templates
 * is that they make the code hard to read and practically all code must be in
 * headers which slows down compilation time.
 *
 * Another approach are virtual methods. Here, one would need to provide functions
 * to create instances of e.g. the Species object. Then, every time the object is
 * retrieved from the library, it needs to be cast to the right class.
 * Using vertial classes means that we must allocate each class and use pointers,
 * this might lead to "scattered?" memory usage compared to vectors.
 *
 * The last approach is to use preprocessor. This brings disadvantages that a big
 * part of code might not be compilable and different implementations can diverge.
 * The code must be maintained in a way that this does not happen and rules must be
 * set.
 *
 * Preprocessor conditioning must be used only in necessary cases such as when
 * defining attributes and methods for Species or Reactions.
 * Everywhere else, configuration from a config class must be used. (TODO: not implemented yet)
 *
 */
#define BNG_ONLY_MCELL


#endif /* LIBS_BNG_CONFIG_H_ */
