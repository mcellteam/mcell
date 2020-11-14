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

#ifndef API_GEN_MODEL_H
#define API_GEN_MODEL_H

#include "../api/common.h"
#include "../api/config.h"
#include "../api/mol_wall_hit_info.h"
#include "../api/notifications.h"
#include "../api/warnings.h"

namespace MCell {
namespace API {

class Config;
class Count;
class ElementaryMoleculeType;
class GeometryObject;
class InstantiationData;
class MolWallHitInfo;
class Molecule;
class Notifications;
class Observables;
class ReactionRule;
class Region;
class ReleaseSite;
class Species;
class Subsystem;
class SurfaceClass;
class VizOutput;
class Wall;
class WallWallHitInfo;
class Warnings;

class GenModel {
public:
  virtual ~GenModel() {}
  std::string to_str(const std::string ind="") const ;

  // --- attributes ---
  Config config;
  virtual void set_config(const Config& new_config_) {
    config = new_config_;
  }
  virtual const Config& get_config() const {
    return config;
  }

  Warnings warnings;
  virtual void set_warnings(const Warnings& new_warnings_) {
    warnings = new_warnings_;
  }
  virtual const Warnings& get_warnings() const {
    return warnings;
  }

  Notifications notifications;
  virtual void set_notifications(const Notifications& new_notifications_) {
    notifications = new_notifications_;
  }
  virtual const Notifications& get_notifications() const {
    return notifications;
  }

  // --- methods ---
  virtual void initialize() = 0;
  virtual void run_iterations(const float_t iterations) = 0;
  virtual void end_simulation(const bool print_final_report = true) = 0;
  virtual void add_subsystem(std::shared_ptr<Subsystem> subsystem) = 0;
  virtual void add_instantiation_data(std::shared_ptr<InstantiationData> instantiation_data) = 0;
  virtual void add_observables(std::shared_ptr<Observables> observables) = 0;
  virtual void dump_internal_state() = 0;
  virtual void export_data_model(const std::string& file = STR_UNSET) = 0;
  virtual void export_viz_data_model(const std::string& file = STR_UNSET) = 0;
  virtual std::vector<int> get_molecule_ids(std::shared_ptr<Species> species = nullptr) = 0;
  virtual std::shared_ptr<Molecule> get_molecule(const int id) = 0;
  virtual Vec3 get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) = 0;
  virtual std::shared_ptr<Wall> get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) = 0;
  virtual Vec3 get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) = 0;
  virtual Vec3 get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) = 0;
  virtual void add_vertex_move(std::shared_ptr<GeometryObject> object, const int vertex_index, const Vec3& displacement) = 0;
  virtual std::vector<std::shared_ptr<WallWallHitInfo>> apply_vertex_moves(const bool collect_wall_wall_hits = false) = 0;
  virtual void register_mol_wall_hit_callback(const std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)> function, py::object context, std::shared_ptr<GeometryObject> object = nullptr, std::shared_ptr<Species> species = nullptr) = 0;
  virtual void load_bngl(const std::string& file_name, const std::string& observables_files_prefix, std::shared_ptr<Region> default_release_region = nullptr, const std::map<std::string, float_t>& parameter_overrides = std::map<std::string, float_t>()) = 0;
  virtual void export_to_bngl(const std::string& file_name) = 0;
}; // GenModel

class Model;
py::class_<Model> define_pybinding_Model(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MODEL_H
