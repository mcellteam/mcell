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

#ifndef API_GEN_COUNT_H
#define API_GEN_COUNT_H

#include "../api/common.h"
#include "../api/count_term.h"


namespace MCell {
namespace API {

class CountTerm;
class GeometryObject;
class ReactionRule;
class Species;

#define COUNT_CTOR() \
    Count( \
        const std::string& filename_, \
        const std::vector<std::shared_ptr<CountTerm>> count_terms_include_ = std::vector<std::shared_ptr<CountTerm>>(), \
        const std::vector<std::shared_ptr<CountTerm>> count_terms_subtract_ = std::vector<std::shared_ptr<CountTerm>>(), \
        const int every_n_timesteps_ = 1, \
        std::shared_ptr<Species> species_ = nullptr, \
        std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, \
        std::shared_ptr<GeometryObject> enclosed_in_object_ = nullptr \
    )  : GenCount(species_,reaction_rule_,enclosed_in_object_) { \
      class_name = "Count"; \
      filename = filename_; \
      count_terms_include = count_terms_include_; \
      count_terms_subtract = count_terms_subtract_; \
      every_n_timesteps = every_n_timesteps_; \
      species = species_; \
      reaction_rule = reaction_rule_; \
      enclosed_in_object = enclosed_in_object_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenCount: public CountTerm {
public:
  GenCount( 
      std::shared_ptr<Species> species_ = nullptr, 
      std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, 
      std::shared_ptr<GeometryObject> enclosed_in_object_ = nullptr 
  )  : CountTerm(species_,reaction_rule_,enclosed_in_object_)  {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::string filename;
  virtual void set_filename(const std::string& new_filename_) {
    filename = new_filename_;
  }
  virtual const std::string& get_filename() const {
    return filename;
  }

  std::vector<std::shared_ptr<CountTerm>> count_terms_include;
  virtual void set_count_terms_include(const std::vector<std::shared_ptr<CountTerm>> new_count_terms_include_) {
    count_terms_include = new_count_terms_include_;
  }
  virtual std::vector<std::shared_ptr<CountTerm>> get_count_terms_include() const {
    return count_terms_include;
  }

  std::vector<std::shared_ptr<CountTerm>> count_terms_subtract;
  virtual void set_count_terms_subtract(const std::vector<std::shared_ptr<CountTerm>> new_count_terms_subtract_) {
    count_terms_subtract = new_count_terms_subtract_;
  }
  virtual std::vector<std::shared_ptr<CountTerm>> get_count_terms_subtract() const {
    return count_terms_subtract;
  }

  int every_n_timesteps;
  virtual void set_every_n_timesteps(const int new_every_n_timesteps_) {
    every_n_timesteps = new_every_n_timesteps_;
  }
  virtual int get_every_n_timesteps() const {
    return every_n_timesteps;
  }

  // --- methods ---
}; // GenCount

class Count;
py::class_<Count> define_pybinding_Count(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COUNT_H
