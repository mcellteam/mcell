/******************************************************************************
 *
 * Copyright (C) 2021 by
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

#include "bngl_exporter.h"

#include <iostream>
#include <sstream>

#include "world.h"
#include "partition.h"
#include "release_event.h"
#include "mol_or_rxn_count_event.h"
#include "vtk_utils.h"

#include "datamodel_defines.h"

using namespace std;

namespace MCell {

void BNGLExporter::clear_temporaries() {
  nfsim_export = false;
  world = nullptr;
}


// returns empty string if everything went well, nonempty string with error message
std::string BNGLExporter::export_to_bngl(
    World* world_,
    const std::string& file_name,
    const API::BNGSimulationMethod simulation_method) {

  clear_temporaries();
  nfsim_export = simulation_method == API::BNGSimulationMethod::NF;
  world = world_;

  // contains all error messages, best effort approach is used where all parts are attempted to be generated
  string err_msg;

  ofstream out;
  out.open(file_name);
  if (!out.is_open()) {
    return "Could not open output file " + file_name + ".";
  }

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  if (p.get_geometry_objects().size() > 1 && nfsim_export) {
    // fail immediately
    return "BNGL export with NFSim is not supported only for models having 1 geometry object.\n";
  }

  // check if there is a single object
  string single_object_name;
  double single_object_volume = 0;
  double single_object_area = 0;
  if (p.get_geometry_objects().size() == 1) {
    const GeometryObject& obj = p.get_geometry_objects()[0];
    single_object_name = obj.name;

    double volume_internal_units = VtkUtils::get_geometry_object_volume(world, obj);
    if (volume_internal_units == FLT_INVALID) {
      return "Compartment object " + obj.name + " is not watertight and its volume cannot be computed.\n";
    }

    single_object_volume = volume_internal_units * pow(world->config.length_unit, 3);
    single_object_area = Geometry::compute_geometry_object_area(p, obj) * pow(world->config.length_unit, 2);
  }

  stringstream parameters;
  stringstream molecule_types;
  stringstream reaction_rules;
  stringstream compartments;
  parameters << BNG::BEGIN_PARAMETERS << "\n";

  parameters << BNG::IND << "# general parameters\n";
  parameters << BNG::IND << BNG::ITERATIONS << " " << world->total_iterations << "\n";
  parameters << BNG::IND << BNG::MCELL_TIME_STEP << " " << f_to_str(world->config.time_unit) << "\n";

  if (single_object_name == BNG::DEFAULT_COMPARTMENT_NAME) {
    parameters << BNG::IND << BNG::MCELL_DEFAULT_COMPARTMENT_VOLUME << " " << f_to_str(single_object_volume) << "\n";
  }

  err_msg += set_compartment_volumes_and_areas();

  err_msg += world->bng_engine.export_to_bngl(
          parameters, molecule_types, compartments, reaction_rules,
          simulation_method == API::BNGSimulationMethod::NF, single_object_volume, single_object_area);

  // seed species
  stringstream seed_species;
  err_msg += export_releases_to_bngl_seed_species(parameters, seed_species);

  stringstream observables;
  err_msg += export_counts_to_bngl_observables(observables);

  parameters << BNG::END_PARAMETERS << "\n";

  out << parameters.str() << "\n";
  out << molecule_types.str() << "\n";
  out << compartments.str() << "\n";
  out << seed_species.str() << "\n";
  out << observables.str() << "\n";

  out << reaction_rules.str() << "\n";

  generate_simulation_action(out, simulation_method);

  out.close();
  return err_msg;
}


std::string BNGLExporter::set_compartment_volumes_and_areas() {

  std::string err_msg;

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  BNG::BNGData& bng_data = world->bng_engine.get_data();

  // first create a mapping of compartments -> geom. objects
  map<BNG::compartment_id_t, const GeometryObject*> compartment_geom_obj_map;

  for (const GeometryObject& obj: p.get_geometry_objects()) {
    if (obj.name == BNG::DEFAULT_COMPARTMENT_NAME) {
      // skip - not generated as compartment
      continue;
    }

    // get compartment information
    BNG::Compartment* compartment = bng_data.find_compartment(obj.name);
    if (compartment == nullptr) {
      err_msg += "Error: geometry object " + obj.name + " is not a BNGL compartment and cannot be exported.\n";
      continue;
    }
    release_assert(compartment->is_3d);

    compartment_geom_obj_map[compartment->id] = &obj;

    if (compartment->parent_compartment_id != BNG::COMPARTMENT_ID_INVALID) {
      BNG::Compartment& parent_compartment = bng_data.get_compartment(compartment->parent_compartment_id);
      if (!parent_compartment.is_3d) {
        // map also surface compartment
        compartment_geom_obj_map[parent_compartment.id] = &obj;
      }
    }
  }

  // now start from compartments that have no children and gradually compute volumes
  std::vector<BNG::compartment_id_t> sorted_compartment_ids;
  bng_data.get_compartments_sorted_by_parents_first(sorted_compartment_ids);

  for (int i = sorted_compartment_ids.size() - 1; i >= 0; i--) {
    BNG::compartment_id_t comp_id = sorted_compartment_ids[i];
    auto it = compartment_geom_obj_map.find(comp_id);
    release_assert(it != compartment_geom_obj_map.end() && "Internal error during compartment export");
    const GeometryObject* obj = it->second;
    BNG::Compartment& comp = bng_data.get_compartment(comp_id);

    if (comp.is_3d) {
      // get volume of all children - we start from compartments with no children so the
      // values were already computed
      comp.set_volume(0);
      double children_volume = comp.get_volume_including_children(bng_data, true);

      // compute total volume
      double volume_internal_units = VtkUtils::get_geometry_object_volume(world, *obj);
      if (volume_internal_units == FLT_INVALID) {
        return "Compartment object " + obj->name + " is not watertight and its volume cannot be computed.\n";
      }
      // must subtract volume of children
      double volume = volume_internal_units * pow(world->config.length_unit, 3) - children_volume;
      comp.set_volume(volume);
    }
    else {
      double area = Geometry::compute_geometry_object_area(p, *obj) * pow(world->config.length_unit, 2);
      comp.set_area(area);
    }
  }
  return err_msg;
}


// returns true if check passed
static std::string get_explicit_compartment_name(
    const BNG::BNGData& bng_data,
    const std::string& err_suffix,
    const bool vol_release,
    const BNG::compartment_id_t release_compartment_id,
    const BNG::compartment_id_t species_compartment_id,
    std::string& explicit_compartment_name
) {
  std::string err_msg;

  const BNG::Compartment& comp = bng_data.get_compartment(release_compartment_id);
  explicit_compartment_name = comp.name;

  if (species_compartment_id != BNG::COMPARTMENT_ID_NONE) {
    // check match
    if (explicit_compartment_name != bng_data.get_compartment(species_compartment_id).name) {
      return string(vol_release ? "Volume" : "Surface") +
          " molecule's compartment does not match the target object compartment" + err_suffix;
    }
    // no need to generate compartment - species has the correct compartment
    explicit_compartment_name = "";
  }

  if (vol_release) {
    if (!comp.is_3d) {
      return "Volume molecule's compartment is a surface/2D compartment" + err_suffix;
    }
  }
  else {
    if (comp.is_3d) {
      return "Surface molecule's compartment is a volume/3D compartment" + err_suffix;
    }
  }

  return "";
}


// search for the leftmost volume compartment if in a tree composed only from differences and leaves
static std::string get_leftmost_compartment_id_recursively(
    const World* world, const std::string& err_suffix, const RegionExprNode* node, BNG::compartment_id_t& id) {

  if (node->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    const GeometryObject& obj = world->get_geometry_object(node->geometry_object_id);
    if (obj.vol_compartment_id == BNG::COMPARTMENT_ID_NONE) {
      return "Trying to convert a volume molecule release for BNGL export but used object " + obj.name + " has "
          "no volume compartment specified" + err_suffix;
    }
    id = obj.vol_compartment_id;
    return "";
  }
  else if (node->op == RegionExprOperator::DIFFERENCE) {
    return get_leftmost_compartment_id_recursively(world, err_suffix, node->left, id);
  }
  else {
    return "Unsupported volume region operator encountered when converting release for BNGL export" + err_suffix;
  }
}


// search for the all volume compartments if in a tree composed only from differences and leaves
static std::string get_all_compartment_ids_recursively_right_is_leaf(
    const World* world, const std::string& err_suffix, const RegionExprNode* node, std::vector<BNG::compartment_id_t>& ids) {

  if (node->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    const GeometryObject& obj = world->get_geometry_object(node->geometry_object_id);
    if (obj.vol_compartment_id == BNG::COMPARTMENT_ID_NONE) {
      return "Trying to convert a volume molecule release for BNGL export but used object " + obj.name + " has "
          "no volume compartment specified" + err_suffix;
    }
    ids.push_back(obj.vol_compartment_id);
    return "";
  }
  else if (node->op == RegionExprOperator::DIFFERENCE) {
    string msg;
    msg = get_all_compartment_ids_recursively_right_is_leaf(world, err_suffix, node->left, ids);
    if (msg != "") {
      return msg;
    }

    if (node->right->op != RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
      return "Unsupported form of region operator expression encountered when converting release for BNGL export" + err_suffix;
    }
    msg = get_all_compartment_ids_recursively_right_is_leaf(world, err_suffix, node->right, ids);
    return msg;
  }
  else {
    return "Unsupported volume region operator encountered when converting release for BNGL export" + err_suffix;
  }
}


std::string BNGLExporter::export_releases_to_bngl_seed_species(
    std::ostream& parameters, std::ostream& seed_species) const {

  string err_msg;

  seed_species << BNG::BEGIN_SEED_SPECIES << "\n";

  parameters << "\n" << BNG::IND << "# seed species counts\n";

  vector<BaseEvent*> release_events;
  world->scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_RELEASE, release_events);

  for (size_t i = 0; i < release_events.size(); i++) {
    ReleaseEvent* re = dynamic_cast<ReleaseEvent*>(release_events[i]);
    assert(re != nullptr);

    if (re->release_shape == ReleaseShape::INITIAL_SURF_REGION) {
      // TODO: check whether there are initial releases, ignoring it for now because it is not often used
      continue;
    }

    const string & err_suffix = ", error for " + re->release_site_name + ".\n";
    if (re->release_shape != ReleaseShape::REGION) {
      err_msg += "Only region release shapes are currently supported for BNGL export" + err_suffix;
      continue;
    }
    if (re->release_number_method != ReleaseNumberMethod::CONST_NUM) {
      err_msg += "Only constant release number releases are currently supported for BNGL export" + err_suffix;
      continue;
    }
    if (re->event_time != 0) {
      err_msg += "Only releases for time 0 are currently supported for BNGL export" + err_suffix;
      continue;
    }
    if (re->orientation != ORIENTATION_NONE && re->orientation != ORIENTATION_UP) {
      err_msg += "Only releases for volume molecules or surface molecules with orientation UP (') "
          "are currently supported for BNGL export" + err_suffix;
      continue;
    }
    if (re->needs_release_pattern()) {
      err_msg += "Releases with release patterns are not currently supported for BNGL export" + err_suffix;
      continue;
    }

    // create parameter for the count
    string seed_count_name = "seed_count_" + to_string(i);
    parameters << BNG::IND << seed_count_name << " " << to_string(re->release_number) << "\n";

    // determine compartment
    const BNG::Species& species = world->bng_engine.get_all_species().get(re->species_id);
    BNG::compartment_id_t species_compartment_id = species.get_primary_compartment_id();
    const BNG::BNGData& bng_data = world->bng_engine.get_data();

    string explicit_compartment_name = "";
    const string err_mgs_only_compartments = "Only release regions that represent a compartment without its children "
        "are supported for BNGL export" + err_suffix;

    if (re->region_expr.root->op == RegionExprOperator::LEAF_SURFACE_REGION) {
      // surface - must be a 2D release

      const Region& region = world->get_region(re->region_expr.root->region_id);
      if (DMUtils::get_region_name(region.name) != "ALL") {
        return "Compartments that do not span the whole object are not supported yet" + err_suffix;
      }
      const GeometryObject& obj = world->get_geometry_object(region.geometry_object_id);

      if (obj.surf_compartment_id == BNG::COMPARTMENT_ID_NONE) {
        err_msg += "Trying to convert a surface molecule release for BNGL export but the object's surface has "
            "no compartment specified" + err_suffix;
        continue;
      }

      string msg = get_explicit_compartment_name(
          bng_data, err_suffix, false, obj.surf_compartment_id, species_compartment_id,
          explicit_compartment_name);

      if (msg != "") {
        err_msg += msg;
        continue;
      }
    }
    else if (re->region_expr.root->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
      const GeometryObject& obj = world->get_geometry_object(re->region_expr.root->geometry_object_id);

      if (obj.name == BNG::DEFAULT_COMPARTMENT_NAME) {
        explicit_compartment_name = "";
      }
      else {
        if (obj.vol_compartment_id == BNG::COMPARTMENT_ID_NONE) {
          err_msg += "Trying to convert a volume molecule release for BNGL export but the object has "
              "no volume compartment specified" + err_suffix;
          continue;
        }

        string msg = get_explicit_compartment_name(
            bng_data, err_suffix, true, obj.vol_compartment_id, species_compartment_id,
            explicit_compartment_name);

        if (msg != "") {
          err_msg += msg;
          continue;
        }
      }
    }
    else if (re->region_expr.root->op == RegionExprOperator::DIFFERENCE) {
      // this may be a volume compartment that has children -> we must check that the
      // release is exactly for this compartment
      // example:
      // EC with CP1 & CP2 -> the region expression must be in this form
      // ((EC - CP1) - CP2)

      BNG::compartment_id_t top_compartment_id;
      string msg = get_leftmost_compartment_id_recursively(world, err_suffix, re->region_expr.root, top_compartment_id);
      if (msg != "") {
        err_msg += msg;
        continue;
      }

      vector<BNG::compartment_id_t> all_compartments;
      msg = get_all_compartment_ids_recursively_right_is_leaf(world, err_suffix, re->region_expr.root, all_compartments);
      if (msg != "") {
        err_msg += msg;
        continue;
      }

      // now check that the compartment corresponds to
      const BNG::Compartment& top_compartment = bng_data.get_compartment(top_compartment_id);

      // number of children must match
      if (all_compartments.size() != top_compartment.children_compartments.size() + 1) {
        err_msg += err_mgs_only_compartments;
        continue;
      }

      set<BNG::compartment_id_t> all_compartments_set(all_compartments.begin(), all_compartments.end());
      // check that the same children are used
      for (BNG::compartment_id_t id: top_compartment.children_compartments) {
        if (all_compartments_set.count(id) != 0) {
          err_msg += err_mgs_only_compartments;
          continue;
        }
      }

      msg = get_explicit_compartment_name(
          bng_data, err_suffix, true, top_compartment_id, species_compartment_id,
          explicit_compartment_name);

      if (msg != "") {
        err_msg += msg;
        continue;
      }
    }
    else {
      err_msg += err_mgs_only_compartments;
      continue;
    }

    if (explicit_compartment_name != "") {
      seed_species << BNG::IND << "@" << explicit_compartment_name << ":";
    }
    else {
      seed_species << BNG::IND;
    }

    seed_species << species.name << " " << seed_count_name << "\n";
  }

  seed_species << BNG::END_SEED_SPECIES << "\n";
  return err_msg;
}


std::string BNGLExporter::export_counts_to_bngl_observables(std::ostream& observables) const {

  observables << BNG::BEGIN_OBSERVABLES << "\n";

  vector<BaseEvent*> count_events;
  world->scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT, count_events);

  for (size_t i = 0; i < count_events.size(); i++) {
    MolOrRxnCountEvent* ce = dynamic_cast<MolOrRxnCountEvent*>(count_events[i]);
    assert(ce != nullptr);

    for (const MolOrRxnCountItem& item: ce->mol_rxn_count_items) {

      const CountBuffer& buff = world->get_count_buffer(item.buffer_id);
      string name;

      // TODO: unify - we should use just column name
      if (buff.get_output_format() == CountOutputFormat::DAT) {
        // get observable name from filename
        const string& path = world->get_count_buffer(item.buffer_id).get_filename();
        size_t slash_pos = path.find_last_of("/\\");
        release_assert(slash_pos != string::npos);
        size_t dot_pos = path.rfind('.');
        release_assert(dot_pos != string::npos);
        name = path.substr(slash_pos + 1, dot_pos - slash_pos - 1);
      }
      else {
        // use column name
        name = buff.get_column_name(item.buffer_column_index);
      }

      const string& err_suffix = ", error for " + name + ".";

      if (item.multiplier != 1.0) {
        return "Observable expressions with a multiplier are not supported by BNGL export yet" + err_suffix;
      }

      string pattern;
      string type;
      for (const MolOrRxnCountTerm& term: item.terms) {
        if (term.sign_in_expression != 1) {
          return "Observable counts with negative value (subtracted) not supported by BNGL export" + err_suffix;
        }
        if (term.is_rxn_count()) {
          return "Reaction counts are not supported by BNGL export" + err_suffix;
        }
        if (term.type == CountType::PresentOnSurfaceRegion) {
          return "Surface molecule counts are not supported by BNGL export yet" + err_suffix;
        }
        if (term.region_expr.root != nullptr && term.region_expr.root->has_binary_op()) {
          return "Counts that use region expressions with union, difference, or intersection are "
              "not supported by BNGL export yet" + err_suffix;
        }

        string compartment_prefix = "";
        if (term.type == CountType::EnclosedInVolumeRegion) {
          release_assert(term.region_expr.root->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT);
          const string& compartment_name = world->get_geometry_object(term.region_expr.root->geometry_object_id).name;
          if (compartment_name != BNG::DEFAULT_COMPARTMENT_NAME) {
            compartment_prefix = "@" + world->get_geometry_object(term.region_expr.root->geometry_object_id).name + ":";
          }
        }

        // Species or Molecules
        string term_type;
        if (term.species_pattern_type == SpeciesPatternType::SpeciesPattern) {
          term_type = BNG::OBSERVABLE_SPECIES;
        }
        else if (term.species_pattern_type == SpeciesPatternType::MoleculesPattern) {
          term_type = BNG::OBSERVABLE_MOLECULES;
        }
        else {
          release_assert("SpeciesId type should not be used here.");
        }
        if (type != "" && type != term_type) {
          return "Combined Molecules and Species observables in one count are not supported by BNGL export" + err_suffix;
        }
        else {
          type = term_type;
        }

        pattern += compartment_prefix + term.species_molecules_pattern.to_str(false, false) + " ";
      }

      observables << BNG::IND <<
          type << " " << name << " " << pattern << " " << "\n";
    }
  }

  observables << BNG::END_OBSERVABLES << "\n";
  return "";
}


void BNGLExporter::generate_simulation_action(
    std::ostream& out, const API::BNGSimulationMethod simulation_method) const {

  if (simulation_method != API::BNGSimulationMethod::NF) {
    out << "generate_network({overwrite=>1})\n";
  }

  string method;
  switch (simulation_method) {
    case API::BNGSimulationMethod::NONE:
      // nothing to do
      break;
    case API::BNGSimulationMethod::ODE:
      method = "ode";
      break;
    case API::BNGSimulationMethod::PLA:
      method = "pla";
      break;
    case API::BNGSimulationMethod::SSA:
      method = "ssa";
      break;
    case API::BNGSimulationMethod::NF:
      method = "nf";
      break;
    default:
      release_assert(false && "Invalid BNG simulation method.");
  }

  if (simulation_method != API::BNGSimulationMethod::NONE) {
    // get sampling frequency from observables
    std::vector<BaseEvent*> count_events;
    world->scheduler.get_all_events_with_type_index(EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT, count_events);
    double min_periodicity = FLT_INVALID;
    for (const BaseEvent* e: count_events) {
      if (e->periodicity_interval != 0 && e->periodicity_interval < min_periodicity) {
        min_periodicity = e->periodicity_interval;
      }
    }
    if (min_periodicity == FLT_INVALID) {
      min_periodicity = 1;
    }

    out <<
        "simulate({method=>\"" << method << "\"," <<
        "seed=>1," <<
        "t_end=>" << world->total_iterations * world->config.time_unit << ","
        "n_steps=>" << world->total_iterations / min_periodicity;

    if (simulation_method == API::BNGSimulationMethod::NF) {
      out << ",glm=>1000000"; // just some default max. molecule count
    }

    out << "})\n";
  }
}

} // namespace MCell
