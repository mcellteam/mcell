/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_GEN_MODEL_H
#define API_GEN_MODEL_H

#include "api/api_common.h"
#include "api/api_config.h"
#include "api/mol_wall_hit_info.h"
#include "api/notifications.h"
#include "api/reaction_info.h"
#include "api/warnings.h"
#include "api/subsystem.h"
#include "api/instantiation.h"
#include "api/observables.h"
#include "api/introspection.h"

namespace MCell {
namespace API {

class Model;
class BaseChkptMol;
class Color;
class Complex;
class Config;
class Count;
class ElementaryMoleculeType;
class GeometryObject;
class Instantiation;
class MolWallHitInfo;
class Molecule;
class Notifications;
class Observables;
class ReactionInfo;
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
class PythonExportContext;

class GenModel: public Subsystem, public Instantiation, public Observables, public Introspection {
public:
  GenModel() {
  }
  GenModel(DefaultCtorArgType) {
  }
  virtual ~GenModel() {}
  std::shared_ptr<Model> copy_model() const;
  std::shared_ptr<Model> deepcopy_model(py::dict = py::dict()) const;
  virtual bool __eq__(const Model& other) const;
  virtual bool eq_nonarray_attributes(const Model& other, const bool ignore_name = false) const;
  bool operator == (const Model& other) const { return __eq__(other);}
  bool operator != (const Model& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const ;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);
  virtual std::string export_vec_species(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_reaction_rules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_surface_classes(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_elementary_molecule_types(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_release_sites(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_geometry_objects(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_checkpointed_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_viz_outputs(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_counts(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


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
  virtual void initialize(const bool print_copyright = true) = 0;
  virtual uint64_t run_iterations(const double iterations) = 0;
  virtual void end_simulation(const bool print_final_report = true) = 0;
  virtual void add_subsystem(std::shared_ptr<Subsystem> subsystem) = 0;
  virtual void add_instantiation(std::shared_ptr<Instantiation> instantiation) = 0;
  virtual void add_observables(std::shared_ptr<Observables> observables) = 0;
  virtual void dump_internal_state(const bool with_geometry = false) = 0;
  virtual void export_data_model(const std::string& file = STR_UNSET) = 0;
  virtual void export_viz_data_model(const std::string& file = STR_UNSET) = 0;
  virtual void export_geometry(const std::string& output_files_prefix = STR_UNSET) = 0;
  virtual void release_molecules(std::shared_ptr<ReleaseSite> release_site) = 0;
  virtual std::vector<int> run_reaction(std::shared_ptr<ReactionRule> reaction_rule, const std::vector<int> reactant_ids, const double time) = 0;
  virtual void add_vertex_move(std::shared_ptr<GeometryObject> object, const int vertex_index, const std::vector<double> displacement) = 0;
  virtual std::vector<std::shared_ptr<WallWallHitInfo>> apply_vertex_moves(const bool collect_wall_wall_hits = false) = 0;
  virtual void register_mol_wall_hit_callback(const std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)> function, py::object context, std::shared_ptr<GeometryObject> object = nullptr, std::shared_ptr<Species> species = nullptr) = 0;
  virtual void register_reaction_callback(const std::function<void(std::shared_ptr<ReactionInfo>, py::object)> function, py::object context, std::shared_ptr<ReactionRule> reaction_rule) = 0;
  virtual void load_bngl(const std::string& file_name, const std::string& observables_path_or_file = STR_UNSET, std::shared_ptr<Region> default_release_region = nullptr, const std::map<std::string, double>& parameter_overrides = std::map<std::string, double>(), const CountOutputFormat observables_output_format = CountOutputFormat::AUTOMATIC_FROM_EXTENSION) = 0;
  virtual void export_to_bngl(const std::string& file_name, const BNGSimulationMethod simulation_method = BNGSimulationMethod::ODE) = 0;
  virtual void save_checkpoint(const std::string& custom_dir = STR_UNSET) = 0;
  virtual void schedule_checkpoint(const uint64_t iteration = 0, const bool continue_simulation = false, const std::string& custom_dir = STR_UNSET) = 0;
}; // GenModel

class Model;
py::class_<Model> define_pybinding_Model(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MODEL_H
