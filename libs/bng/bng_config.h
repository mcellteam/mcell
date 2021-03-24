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

// separate class because the number of notifications may grow in the future
class BNGNotifications {
public:
  BNGNotifications() :
    bng_verbosity_level(0),
    rxn_probability_changed(true)
    {
  }

  // print additional information during simulation
  int bng_verbosity_level;

  // same as MCell3's varying_probability_report,
  // TODO: probably not converted from data model
  bool rxn_probability_changed;

  void dump() const;
};


class BNGWarnings {
public:
  BNGWarnings() :
    warn_on_bimol_rxn_probability_over_05_less_1(true),
    bimol_rxn_probability_over_05_less_1(false),
    bimol_rxn_probability_over_1(false)
    {
  }

  // warning settings
  bool warn_on_bimol_rxn_probability_over_05_less_1;

  // status
  mutable bool bimol_rxn_probability_over_05_less_1; // set to true if probability of a bimol rxn was > 0.5 but < 1
  mutable bool bimol_rxn_probability_over_1; // set to true if probability of a bimol rxn was > 1

};


class BNGConfig {
public:
  BNGConfig() :
    mcell(BNG_MCELL),
    initial_seed(1),
    time_unit(0),
    length_unit(0),
    rcp_length_unit(0),
    grid_density(0),
    rx_radius_3d(0),
    intermembrane_rx_radius_3d(0),
    rxn_and_species_report(true)
    {
  }

  // configuration
  bool mcell;

  uint initial_seed;

  // MCell
  BNGCommon::float_t time_unit;
  BNGCommon::float_t length_unit;
  BNGCommon::float_t rcp_length_unit;
  BNGCommon::float_t grid_density;
  BNGCommon::float_t rx_radius_3d;
  BNGCommon::float_t intermembrane_rx_radius_3d;

  // generate report files during simulation
  bool rxn_and_species_report;

  BNGWarnings warnings;
  BNGNotifications notifications;


  void init() {
    rcp_length_unit = 1/length_unit;
    initialize_bng_report_files();
  }

  std::string get_rxn_report_file_name() const;
  std::string get_species_report_file_name() const;
  std::string get_warnings_report_file_name() const;

  // prints warnings at the end of simulation if needed (to the standard error output)
  void print_final_warnings() const;

  void dump() const;

  BNGCommon::float_t get_default_rx_radius_3d() const;

protected:
  std::string seed_as_str() const;
  void initialize_bng_report_files();
};


// report utilities
std::string get_current_date_time();
void append_to_report(const std::string& report_fname, const std::string& msg);
void initialize_report_file(const std::string& fname, const char* report_name);
void remove_report_file(const std::string& fname, const char* report_name);

} // namespace BNG


#endif /* LIBS_BNG_BNG_CONFIG_H_ */
