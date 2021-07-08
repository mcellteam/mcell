/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef LIBMCELL_API_PYTHON_EXPORTER_H_
#define LIBMCELL_API_PYTHON_EXPORTER_H_

#include <string>
#include <map>
#include <ostream>

namespace MCell {

class World;

namespace API {

class Model;
class PythonExportContext;

class PythonExporter {
public:
  PythonExporter(Model* model_);

  void save_checkpoint(const std::string& output_dir_);
private:
  void open_and_check_file(
      const std::string file_name, std::ofstream& out,
      const bool for_append = false,
      const bool bngl = false);

  std::string export_subsystem(PythonExportContext& ctx);
  std::string export_geometry(PythonExportContext& ctx);
  std::string export_instantiation(PythonExportContext& ctx, const std::string& geometry_objects_name);
  std::string export_observables(PythonExportContext& ctx);
  void export_simulation_state(
      PythonExportContext& ctx, std::map<std::string, std::string>& config_variable_names);
  void export_molecules(std::ostream& out, PythonExportContext& ctx);

  std::string export_model(
      PythonExportContext& ctx,
      const std::string& subsystem_name,
      const std::string& instantiation_name,
      const std::string& observables_name,
      const std::map<std::string, std::string>& config_variable_names);
  void export_checkpoint_iterations(std::ostream& out);


  Model* model;
  World* world;

  std::string output_dir;
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_PYTHON_EXPORTER_H_ */
