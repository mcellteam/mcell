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

#include <exception>
#include <string>

#include "bng_data_to_datamodel_converter.h"

#include "datamodel_defines.h"
#include "bng/bng.h"

#include "viz_output_event.h"
#include "scheduler.h"
#include "world.h"

using namespace std;
using namespace DMUtils;
using namespace BNG;
using Json::Value;

namespace MCell {

typedef std::invalid_argument ConversionError;

#define CHECK(stmt) \
  do { \
    try { \
      (stmt); \
    } \
    catch (exception& e) { \
      cerr << e.what() << "\n"; \
      cerr << "Exception caught in '" << __FUNCTION__ << "' after conversion error.\n"; \
      conversion_failed = true; \
    } \
  } while (0)

#define CHECK_PROPERTY(cond) \
  do { \
    if(!(cond)) { \
      throw ConversionError("Expected '" + #cond + "' is false. (" + __FUNCTION__ + " - " + __FILE__ + ":" + to_string(__LINE__) + ")"); \
    } \
  } while (0)

#define CHECK_PROPERTY_W_NAME(cond, name) \
  do { \
    if(!(cond)) { \
      throw ConversionError(string("Conversion of ") + name + ": Expected '" + #cond + "' is false. (" + __FUNCTION__ + " - " + __FILE__ + ":" + to_string(__LINE__) + ")"); \
    } \
  } while (0)

BngDataToDatamodelConverter::BngDataToDatamodelConverter() {
  reset();

  // list from https://www.rapidtables.com/web/color/RGB_Color.html
  colors.push_back(Vec3(1, 0, 0));
  colors.push_back(Vec3(0, 1, 0));
  colors.push_back(Vec3(0, 0, 1));
  colors.push_back(Vec3(1, 1, 0));
  colors.push_back(Vec3(0, 1, 1));
  colors.push_back(Vec3(1, 0, 1));
  // ...
}


void BngDataToDatamodelConverter::reset() {
  world = nullptr;
  bng_engine = nullptr;
  next_color_index = 0;
  rxn_counter = 0;
  conversion_failed = false;
  processed_surface_classes.clear();
}


void BngDataToDatamodelConverter::to_data_model(
    const World* world_, Value& mcell_node, const bool only_for_viz) {

  reset();

  world = world_;
  bng_engine = &(world->bng_engine);

  // similarly as in pymcell converter, maximum effort is given to conversion and
  // if some parts won't go through, we are trying to generate
  CHECK(convert_molecules(mcell_node));

  // compartments are converted in GeometryObject into the KEY_MODEL_OBJECTS section

  if (!only_for_viz) {
    CHECK(convert_rxns(mcell_node));
  }

  return;
}


Vec3 BngDataToDatamodelConverter::get_next_color() {
  Vec3 res = colors[next_color_index];
  next_color_index++;
  if (next_color_index == colors.size()) {
    next_color_index = 0;
  }
  return res;
}


void BngDataToDatamodelConverter::convert_molecules(Value& mcell_node) {
  Value& define_molecules = mcell_node[KEY_DEFINE_MOLECULES];
  add_version(define_molecules, VER_DM_2014_10_24_1638);
  Value& molecule_list = define_molecules[KEY_MOLECULE_LIST];
  molecule_list = Value(Json::arrayValue); // empty array

  // we are converting molecule types because that's the information we need in
  // the datamodel, we do not care about species
  for (const ElemMolType& mt: bng_engine->get_data().get_elem_mol_types()) {
    // - ALL_* species are ignored
    // - reactive surfaces are converted as surface classes
    // -
    if (is_species_superclass(mt.name) || mt.is_reactive_surface()) {
      continue;
    }
    Value molecule_node;
    CHECK(convert_single_mol_type(mt, molecule_node));
    molecule_list.append(molecule_node);
  }
}


void BngDataToDatamodelConverter::convert_single_mol_type(const BNG::ElemMolType& mt, Value& molecule_node) {

  add_version(molecule_node, VER_DM_2018_10_16_1632);

  Value& display = molecule_node[KEY_DISPLAY];
  display[KEY_EMIT] = 0.0;
  Value& color_value = display[KEY_COLOR];

  if (mt.color_set) {
    color_value.append(mt.color_r);
    color_value.append(mt.color_g);
    color_value.append(mt.color_b);
  }
  else {
    // automatic color
    Vec3 c = get_next_color();
    color_value.append(c.r);
    color_value.append(c.g);
    color_value.append(c.b);
  }
  display[KEY_GLYPH] = VALUE_GLYPH_SPHERE_1;
  display[KEY_SCALE] = mt.scale;

  molecule_node[KEY_BNGL_COMPONENT_LIST] = Value(Json::arrayValue); // empty array
  molecule_node[KEY_MOL_BNGL_LABEL] = "";
  molecule_node[KEY_DESCRIPTION] = "";
  molecule_node[KEY_MOL_TYPE] = mt.is_vol() ? VALUE_MOL_TYPE_3D : VALUE_MOL_TYPE_2D;

  const VizOutputEvent* viz =
      dynamic_cast<const VizOutputEvent*>(world->scheduler.find_next_event_with_type_index(EVENT_TYPE_INDEX_VIZ_OUTPUT));
  if (viz == nullptr || viz->should_visualize_all_species()) {
    // no visualization/or enabled globally
    molecule_node[KEY_EXPORT_VIZ] = false;
  }
  else {
    // simple species based on this molecule type set to be visualized?
    // the simple species have always the same name
    species_id_t species_id = bng_engine->get_all_species().find_by_name(mt.name);
    molecule_node[KEY_EXPORT_VIZ] = viz->species_ids_to_visualize.count(species_id) == 1;
  }

  if (mt.custom_space_step != 0) {
    molecule_node[KEY_CUSTOM_SPACE_STEP] = f_to_str(mt.custom_space_step);
  }
  else {
    molecule_node[KEY_CUSTOM_SPACE_STEP] = "";
  }

  if (mt.custom_time_step != 0) {
    molecule_node[KEY_CUSTOM_TIME_STEP] = f_to_str(mt.custom_time_step);
  }
  else {
    molecule_node[KEY_CUSTOM_TIME_STEP] = "";
  }

  molecule_node[KEY_MAXIMUM_STEP_LENGTH] = "";
  molecule_node[KEY_TARGET_ONLY] = mt.cant_initiate();
  molecule_node[KEY_DIFFUSION_CONSTANT] = f_to_str(mt.D);
  molecule_node[KEY_SPATIAL_STRUCTURE] = "None";
  molecule_node[KEY_MOL_NAME] = mt.name;

  if (!mt.component_type_ids.empty()) {
    // generate info on components
    Value& bngl_component_list = molecule_node[KEY_BNGL_COMPONENT_LIST];
    bngl_component_list = Value(Json::arrayValue);

    for (component_type_id_t ct_id: mt.component_type_ids) {
      Value bngl_component;

      const ComponentType& ct = bng_engine->get_data().get_component_type(ct_id);
      bngl_component[KEY_CNAME] = ct.name;

      Value& cstates = bngl_component[KEY_CSTATES];
      cstates = Value(Json::arrayValue);
      for (state_id_t s_id: ct.allowed_state_ids) {
        const string& state_name = bng_engine->get_data().get_state_name(s_id);
        cstates.append(state_name);
      }

      // additional empty data required by cellblender
      bngl_component[KEY_IS_KEY] = false;
      bngl_component[KEY_ROT_INDEX] = 0;
      bngl_component[KEY_ROT_ANG] = "0";
      bngl_component[KEY_ROT_X] = "0";
      bngl_component[KEY_ROT_Y] = "0";
      bngl_component[KEY_ROT_Z] = "0";
      bngl_component[KEY_LOC_X] = "0";
      bngl_component[KEY_LOC_Y] = "0";
      bngl_component[KEY_LOC_Z] = "0";

      bngl_component_list.append(bngl_component);
    }
  }
}


void BngDataToDatamodelConverter::convert_single_rxn_rule(const BNG::RxnRule& r, Value& rxn_node) {
  assert(r.type == RxnType::Standard);
  add_version(rxn_node, VER_DM_2018_01_11_1330);

  // name is put into rnx_name and name is the string of the reaction
  rxn_node[KEY_RXN_NAME] = r.name;
  rxn_node[KEY_DESCRIPTION] = "";

  // LATER: maybe find an opposite reaction and generate it as reversible
  rxn_node[KEY_RXN_TYPE] = VALUE_IRREVERSIBLE;

  string reactants = r.reactants_to_str();
  rxn_node[KEY_REACTANTS] = reactants;
  string products = r.products_to_str();
  if (products == "") {
    products = VALUE_NULL;
  }
  rxn_node[KEY_PRODUCTS] = products;
  // generated directly by data_model_to_mdl converter, 
  // reversible rxns were cloned
  rxn_node[KEY_NAME] = reactants + " -> " + products;

  if (r.base_variable_rates.empty()) {
    rxn_node[KEY_FWD_RATE] = f_to_str(r.base_rate_constant);
    rxn_node[KEY_BKWD_RATE] = "";

    CHECK_PROPERTY_W_NAME(r.base_variable_rates.empty(), r.name);
    rxn_node[KEY_VARIABLE_RATE_VALID] = false;
    rxn_node[KEY_VARIABLE_RATE_SWITCH] = false;
    rxn_node[KEY_VARIABLE_RATE] = "";
    rxn_node[KEY_VARIABLE_RATE_TEXT] = "";
  }
  else {
    rxn_node[KEY_FWD_RATE] = "0";
    rxn_node[KEY_BKWD_RATE] = "";

    rxn_node[KEY_VARIABLE_RATE_VALID] = true;
    rxn_node[KEY_VARIABLE_RATE_SWITCH] = true;

    rxn_node[KEY_VARIABLE_RATE] = "var_rate_" + DMUtils::to_id(r.name) + "_" + to_string(rxn_counter);
    rxn_counter++;

    stringstream text;
    // initial value for time 0
    text << 0.0 << "\t" << r.base_rate_constant << "\n";
    for (const RxnRateInfo& ri: r.base_variable_rates) {
      text << ri.time * world->config.time_unit  << "\t" << ri.rate_constant << "\n";
    }
    rxn_node[KEY_VARIABLE_RATE_TEXT] = text.str();
  }
}


string BngDataToDatamodelConverter::get_surface_class_name(const BNG::RxnRule& r) {
  assert(r.is_bimol() && r.is_reactive_surface_rxn());
  
  uint idx = (r.reactants[0].is_reactive_surface()) ? 0 : 1;
  assert(r.reactants[idx].is_simple());
  // the name of the surface class is the name of the species/elem mol type
  return
      bng_engine->get_data().get_elem_mol_type(r.reactants[idx].get_simple_species_mol_type_id()).name;
}


void BngDataToDatamodelConverter::convert_single_surface_class(const BNG::RxnRule& base_rxn, Value& surface_class) {
  surface_class[KEY_DESCRIPTION] = "";
  add_version(surface_class, VER_DM_2018_01_11_1330);

  string name = get_surface_class_name(base_rxn);
  assert(processed_surface_classes.count(name) == 0 && "This surface class was already processed");
  processed_surface_classes.insert(name);
  surface_class[KEY_NAME] = name;

  Value& surface_class_prop_list = surface_class[KEY_SURFACE_CLASS_PROP_LIST];
  surface_class_prop_list = Value(Json::arrayValue);

  // find all rxn rules with the same name
  vector<const RxnRule*> rxns_belonging_to_surf_class;
  for (const RxnRule* rxn_rule: bng_engine->get_all_rxns().get_rxn_rules_vector()) {
    if (rxn_rule->is_reactive_surface_rxn() && get_surface_class_name(*rxn_rule) == name) {
      rxns_belonging_to_surf_class.push_back(rxn_rule);
    }
  }

  // and convert them all as properties of this surface class
  for (const RxnRule* rxn_rule: rxns_belonging_to_surf_class) {
    Value sc_property;
    add_version(sc_property, VER_DM_2015_11_08_1756);
    CONVERSION_CHECK(rxn_rule->reactants[0].is_simple(), "Surface class reactant must be simple for now.");

    bool has_type = true;
    switch (rxn_rule->type) {
      case RxnType::Transparent:
        sc_property[KEY_SURF_CLASS_TYPE] = VALUE_TRANSPARENT;
        break;
      case RxnType::Reflect:
        sc_property[KEY_SURF_CLASS_TYPE] = VALUE_REFLECTIVE;
        break;
      case RxnType::AbsorbRegionBorder:
        sc_property[KEY_SURF_CLASS_TYPE] = VALUE_ABSORPTIVE;
        break;
      case RxnType::Standard:
        if (rxn_rule->is_absorptive_region_rxn()) {
          sc_property[KEY_SURF_CLASS_TYPE] = VALUE_ABSORPTIVE;
        }
        else {
          // do not generate any property if type is nto set
          continue;
        }
        break;
      default:
        CONVERSION_UNSUPPORTED("Unexpected reaction type");
    }

    const string& reactant_name = bng_engine->get_data().get_elem_mol_type(rxn_rule->reactants[0].get_simple_species_mol_type_id()).name;
    if (is_species_superclass(reactant_name)) {
      sc_property[KEY_AFFECTED_MOLS] = reactant_name;
      sc_property[KEY_MOLECULE] = "";
    }
    else {
      sc_property[KEY_AFFECTED_MOLS] = VALUE_SINGLE;
      sc_property[KEY_MOLECULE] = reactant_name;
    }

    sc_property[KEY_SURF_CLASS_ORIENT] = DMUtils::orientation_to_str(rxn_rule->reactants[0].get_orientation());

    sc_property[KEY_CLAMP_VALUE] = "0";
    sc_property[KEY_NAME] = ""; // name is ignored by the datamodel to mdl converter anyway

    surface_class_prop_list.append(sc_property);
  }
}


void BngDataToDatamodelConverter::convert_rxns(Value& mcell_node) {

  // --- define_reactions --- (this section is always present)
  Value& define_reactions = mcell_node[KEY_DEFINE_REACTIONS];
  add_version(define_reactions, VER_DM_2014_10_24_1638);
  Value& reaction_list = define_reactions[KEY_REACTION_LIST];
  reaction_list = Value(Json::arrayValue);

  // --- define_surface_classes --- (only when needed)
  Value& define_surface_classes = mcell_node[KEY_DEFINE_SURFACE_CLASSES];
  add_version(define_surface_classes, VER_DM_2014_10_24_1638);
  Value& surface_class_list = define_surface_classes[KEY_SURFACE_CLASS_LIST];
  if (!surface_class_list.isArray()) {
    // do not overwrite if it was already set
    surface_class_list = Value(Json::arrayValue);
  }

  for (const RxnRule* rxn_rule: bng_engine->get_all_rxns().get_rxn_rules_vector()) {
    if (rxn_rule->has_flag(RXN_FLAG_CREATED_FOR_CONCENTRATION_CLAMP) ||
        rxn_rule->has_flag(RXN_FLAG_CREATED_FOR_FLUX_CLAMP) ) {
      // concentration clamp info is generated from ConcentrationClampReleaseEvent
      continue;
    }

    if (rxn_rule->type == RxnType::Standard) {
      // this might be a special case of absorptive reaction
      if (!rxn_rule->is_absorptive_region_rxn()) {
        Value rxn_node;
        CHECK(convert_single_rxn_rule(*rxn_rule, rxn_node));
        reaction_list.append(rxn_node);
      }

      // also register the surface class
      if (rxn_rule->is_reactive_surface_rxn() &&
          processed_surface_classes.count(get_surface_class_name(*rxn_rule)) == 0) {
        Value surface_class;
        CHECK(convert_single_surface_class(*rxn_rule, surface_class));
        surface_class_list.append(surface_class);
      }
    }
    else if (rxn_rule->type == RxnType::Transparent ||
        rxn_rule->type == RxnType::Reflect ||
        rxn_rule->type == RxnType::AbsorbRegionBorder) {
      // this rxn rule defines a surface class
      if (processed_surface_classes.count(get_surface_class_name(*rxn_rule)) == 0) {
        Value surface_class;
        CHECK(convert_single_surface_class(*rxn_rule, surface_class));
        surface_class_list.append(surface_class);
      }
    }
    else {
      CONVERSION_UNSUPPORTED("Invalid reaction type");
    }
  }
}


} // namespace MCell

