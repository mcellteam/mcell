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

#ifndef API_GEN_INSTANTIATION_H
#define API_GEN_INSTANTIATION_H

#include "api/api_common.h"
#include "api/base_export_class.h"

namespace MCell {
namespace API {

class Instantiation;
class BaseChkptMol;
class GeometryObject;
class Region;
class ReleaseSite;
class PythonExportContext;

#define INSTANTIATION_CTOR() \
    Instantiation( \
        const std::vector<std::shared_ptr<ReleaseSite>> release_sites_ = std::vector<std::shared_ptr<ReleaseSite>>(), \
        const std::vector<std::shared_ptr<GeometryObject>> geometry_objects_ = std::vector<std::shared_ptr<GeometryObject>>(), \
        const std::vector<std::shared_ptr<BaseChkptMol>> checkpointed_molecules_ = std::vector<std::shared_ptr<BaseChkptMol>>() \
    ) { \
      release_sites = release_sites_; \
      geometry_objects = geometry_objects_; \
      checkpointed_molecules = checkpointed_molecules_; \
    } \
    Instantiation(DefaultCtorArgType){ \
    }

class GenInstantiation: public BaseExportClass {
public:
  GenInstantiation() {
  }
  GenInstantiation(DefaultCtorArgType) {
  }
  virtual ~GenInstantiation() {}
  std::shared_ptr<Instantiation> copy_instantiation() const;
  std::shared_ptr<Instantiation> deepcopy_instantiation(py::dict = py::dict()) const;
  virtual bool __eq__(const Instantiation& other) const;
  virtual bool eq_nonarray_attributes(const Instantiation& other, const bool ignore_name = false) const;
  bool operator == (const Instantiation& other) const { return __eq__(other);}
  bool operator != (const Instantiation& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const ;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);
  virtual std::string export_vec_release_sites(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_geometry_objects(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_checkpointed_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::shared_ptr<ReleaseSite>> release_sites;
  virtual void set_release_sites(const std::vector<std::shared_ptr<ReleaseSite>> new_release_sites_) {
    release_sites = new_release_sites_;
  }
  virtual std::vector<std::shared_ptr<ReleaseSite>>& get_release_sites() {
    return release_sites;
  }

  std::vector<std::shared_ptr<GeometryObject>> geometry_objects;
  virtual void set_geometry_objects(const std::vector<std::shared_ptr<GeometryObject>> new_geometry_objects_) {
    geometry_objects = new_geometry_objects_;
  }
  virtual std::vector<std::shared_ptr<GeometryObject>>& get_geometry_objects() {
    return geometry_objects;
  }

  std::vector<std::shared_ptr<BaseChkptMol>> checkpointed_molecules;
  virtual void set_checkpointed_molecules(const std::vector<std::shared_ptr<BaseChkptMol>> new_checkpointed_molecules_) {
    checkpointed_molecules = new_checkpointed_molecules_;
  }
  virtual std::vector<std::shared_ptr<BaseChkptMol>>& get_checkpointed_molecules() {
    return checkpointed_molecules;
  }

  // --- methods ---
  virtual void add_release_site(std::shared_ptr<ReleaseSite> s) = 0;
  virtual std::shared_ptr<ReleaseSite> find_release_site(const std::string& name) = 0;
  virtual void add_geometry_object(std::shared_ptr<GeometryObject> o) = 0;
  virtual std::shared_ptr<GeometryObject> find_geometry_object(const std::string& name) = 0;
  virtual std::shared_ptr<GeometryObject> find_volume_compartment_object(const std::string& name) = 0;
  virtual std::shared_ptr<GeometryObject> find_surface_compartment_object(const std::string& name) = 0;
  virtual void load_bngl_compartments_and_seed_species(const std::string& file_name, std::shared_ptr<Region> default_release_region = nullptr, const std::map<std::string, double>& parameter_overrides = std::map<std::string, double>()) = 0;
}; // GenInstantiation

class Instantiation;
py::class_<Instantiation> define_pybinding_Instantiation(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_INSTANTIATION_H
