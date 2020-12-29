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

#ifndef LIBMCELL_API_PYTHON_EXPORTER_H_
#define LIBMCELL_API_PYTHON_EXPORTER_H_

#include <string>
#include <map>

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

  std::string save_subsystem(PythonExportContext& ctx);
  std::string save_geometry(PythonExportContext& ctx);
  std::string save_instantiation(PythonExportContext& ctx, const std::string& geometry_objects_name);
  std::string save_observables(PythonExportContext& ctx);
  void save_simulation_state(
      PythonExportContext& ctx, std::map<std::string, std::string>& config_variable_names);
  void save_molecules(PythonExportContext& ctx, std::ostream& out);
  std::string save_model(
      PythonExportContext& ctx,
      const std::string& subsystem_name,
      const std::string& instantiation_name,
      const std::string& observables_name,
      const std::map<std::string, std::string>& config_variable_names);


  Model* model;
  World* world;

  std::string output_dir;
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_PYTHON_EXPORTER_H_ */
