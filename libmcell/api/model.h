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
#include "api/instantiation_data.h"
#include "api/observables.h"
#include "api/config.h"
#include "api/warnings.h"
#include "api/notifications.h"
#include "api/shared_structs.h"
#include "api/callbacks.h"

namespace MCell {

class World;

namespace API {

class MolWallHitInfo;
class WallWallHitInfo;

class Model: public GenModel {
public:

  Model() : initialized(false), world(nullptr), callbacks(this) {
    // add species superclasses
    species.push_back(AllMolecules);
    species.push_back(AllVolumeMolecules);
    species.push_back(AllSurfaceMolecules);
  }
  virtual ~Model();

  // from generated template
  void initialize() override;
  void run_iterations(const float_t iterations) override;
  void end_simulation(const bool print_final_report = true) override;

  void add_subsystem(std::shared_ptr<Subsystem> subsystem) override;
  void add_instantiation_data(std::shared_ptr<InstantiationData> instantiation_data) override;
  void add_observables(std::shared_ptr<Observables> observables) override;

  void dump_internal_state() override;

  void export_data_model(const std::string& file = STR_UNSET) override {
    export_data_model_viz_or_full(file, false, NAME_EXPORT_DATA_MODEL);
  }
  void export_viz_data_model(const std::string& file = STR_UNSET) override {
    export_data_model_viz_or_full(file, true, NAME_EXPORT_VIZ_DATA_MODEL);
  }

  void release_molecules(std::shared_ptr<ReleaseSite> release_site) override;

  std::vector<int> get_molecule_ids(std::shared_ptr<Species> species = nullptr) override;
  std::shared_ptr<Molecule> get_molecule(const int id) override;

  Vec3 get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) override;
  std::shared_ptr<Wall> get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) override;
  Vec3 get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) override;
  Vec3 get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) override;


  void add_vertex_move(
      std::shared_ptr<GeometryObject> object, const int vertex_index, const Vec3& displacement
  ) override;

  std::vector<std::shared_ptr<WallWallHitInfo>> apply_vertex_moves(const bool collect_wall_wall_hits = false) override;

  void register_mol_wall_hit_callback(
      const std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)> function,
      py::object context,
      std::shared_ptr<GeometryObject> object = nullptr,
      std::shared_ptr<Species> species = nullptr
  ) override;

  void register_reaction_callback(
      const std::function<void(std::shared_ptr<ReactionInfo>, py::object)> function,
      py::object context,
      std::shared_ptr<ReactionRule> reaction_rule
  ) override;

  void load_bngl(
      const std::string& file_name,
      const std::string& observables_files_prefix = "",
      std::shared_ptr<Region> default_release_region = nullptr,
      const std::map<std::string, float_t>& parameter_overrides = std::map<std::string, float_t>()
  ) override;

  void export_to_bngl(const std::string& file_name) override;

  // overrides from derived classes Subsystem, InstantiationData, and Observables
  void add_species(std::shared_ptr<Species> s) override;
  void add_reaction_rule(std::shared_ptr<ReactionRule> r) override;

  void add_surface_class(std::shared_ptr<SurfaceClass> sc) override;
  void add_release_site(std::shared_ptr<ReleaseSite> s) override;
  void add_geometry_object(std::shared_ptr<GeometryObject> o) override;

  void add_viz_output(std::shared_ptr<VizOutput> viz_output) override;
  void add_count(std::shared_ptr<Count> count) override;

  // added manually
  // shadows all inherited non-virtual to_str methods
  std::string to_str(const std::string ind="") const;

  std::shared_ptr<GeometryObject> get_geometry_object_with_id(const geometry_object_id_t id);
  std::shared_ptr<ReactionRule> get_reaction_rule_with_fwd_id(const BNG::rxn_rule_id_t id);

  void error_if_initialized(const char* what) {
    if (initialized) {
      throw RuntimeError(S("It is not possible to add ") + what + " once a model was initialized.");
    }
  }

  void dump() const;

  const World* get_world() const {
    return world;
  }

private:
  void export_data_model_viz_or_full(
      const std::string& file,
      const bool only_for_visualization,
      const char* method_name
  );

  bool initialized;
  World* world;

  Callbacks callbacks;

  std::vector<VertexMoveInfo> vertex_moves;
};

} // namespace API
} // namespace MCell

#endif // API_MODEL_H
