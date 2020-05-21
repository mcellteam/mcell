#ifndef LIBS_BNG_BNG_CONFIG_H_
#define LIBS_BNG_BNG_CONFIG_H_

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
 * One more approach is to use preprocessor. This brings disadvantages that a big
 * part of code might not be compilable and different implementations can diverge.
 * For example one could use preprocessor conditioning must be used only in necessary
 * cases such as when defining attributes and methods for Species or Reactions,
 * however,
 *
 * Therefore, all the variables needed for species and reactions along with their
 * code for every simulator that uses this library is present here and shared.
 * We might rethink this later.
 *
 */
#define BNG_MCELL 1

// included only in defines_shared.h

namespace BNG {

class BNGNotifications {
public:
  BNGNotifications()
    : rxn_probability_changed(true)
    {
  }

  bool rxn_probability_changed; // related to MCell's varying_probability_report?

  void dump() const;
};

class BNGConfig {
public:
  BNGConfig()
    : mcell(BNG_MCELL),

      time_unit(0),
      length_unit(0),
      grid_density(0),
      rx_radius_3d(0),

      debug_reactions(false)
  {
  }

  // configuration
  bool mcell;

  // MCell
  BNGCommon::float_t time_unit;
  BNGCommon::float_t length_unit;
  BNGCommon::float_t grid_density;
  BNGCommon::float_t rx_radius_3d;

  // debug
  bool debug_reactions;

  BNGNotifications notifications;

  void dump() const;
};

} // namespace BNG


#endif /* LIBS_BNG_BNG_CONFIG_H_ */