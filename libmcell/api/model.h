/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#ifndef API_MODEL_H
#define API_MODEL_H

#include "generated/gen_model.h"
#include "api/common.h"
#include "api/globals.h"
#include "api/subsystem.h"
#include "api/instantiation_data.h"
#include "api/observables.h"
#include "api/config.h"
#include "api/warnings.h"
#include "api/notifications.h"
#include "api/shared_structs.h"

namespace MCell {

class World;

namespace API {

class Model: public GenModel, public Subsystem, public InstantiationData, public Observables {
public:

  Model() : initialized(false), world(nullptr) {
    // add species superclasses
    species.push_back(AllMolecules);
    species.push_back(AllVolumeMolecules);
    species.push_back(AllSurfaceMolecules);
  }
  virtual ~Model();

  // from generated template
  void initialize() override;
  void run_iterations(const long iterations) override;
  void end_simulation(const bool print_final_report = true) override;

  void add_subsystem(std::shared_ptr<Subsystem> subsystem) override;
  void add_instantiation_data(std::shared_ptr<InstantiationData> instantiation_data) override;
  void add_observables(std::shared_ptr<Observables> observables) override;

  void dump_internal_state() override;

  void export_data_model(const std::string& file = STR_UNSET, const bool only_for_visualization = false) override;

  std::vector<int> get_molecule_ids(std::shared_ptr<Species> species = nullptr) override;

  std::shared_ptr<Molecule> get_molecule(const int id) override;

  void add_vertex_move(
      std::shared_ptr<GeometryObject> object, const int index, const Vec3& displacement
  ) override;

  void apply_vertex_moves() override;

  void register_wall_hit_callback(
      const std::function<void(std::shared_ptr<WallHitInfo>, py::object)> function,
      py::object context,
      std::shared_ptr<GeometryObject> object = nullptr,
      std::shared_ptr<Species> species = nullptr
  ) override;

  void error_if_initialized(const char* what) {
    if (initialized) {
      throw(S("It is not possible to add ") + what + " once a model was initialized.");
    }
  }

  // overrides from derived classes Subsystem, InstantiationData, and Observables
  void add_species(std::shared_ptr<Species> s) override {
    error_if_initialized(NAME_CLASS_SPECIES);
    Subsystem::add_species(s);
  }

  void add_reaction_rule(std::shared_ptr<ReactionRule> r) override {
    error_if_initialized(NAME_CLASS_REACTION_RULE);
    Subsystem::add_reaction_rule(r);
  }

  void add_surface_class(std::shared_ptr<SurfaceClass> sc) override {
    error_if_initialized(NAME_CLASS_SURFACE_CLASS);
    Subsystem::add_surface_class(sc);
  }

  void add_release_site(std::shared_ptr<ReleaseSite> s) override {
    error_if_initialized(NAME_CLASS_RELEASE_SITE);
    InstantiationData::add_release_site(s);
  }

  void add_geometry_object(std::shared_ptr<GeometryObject> o) override {
    error_if_initialized(NAME_CLASS_GEOMETRY_OBJECT);
    InstantiationData::add_geometry_object(o);
  }

  void add_viz_output(std::shared_ptr<VizOutput> viz_output) override {
    error_if_initialized(NAME_CLASS_VIZ_OUTPUT);
    Observables::add_viz_output(viz_output);
  };

  void add_count(std::shared_ptr<Count> count) override {
    error_if_initialized(NAME_CLASS_OBSERVABLES);
    Observables::add_count(count);
  };

  // added manually
  // shadows all inherited non-virtual to_str methods
  std::string to_str(const std::string ind="") const;

  void dump() const;

private:
  bool initialized;
  World* world;

  std::vector<VertexMoveInfo> vertex_moves;
};

} // namespace API
} // namespace MCell

#endif // API_MODEL_H
