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

#include <exception>
#include <string>

#include "bng_data_to_datamodel_converter.h"

#include "datamodel_defines.h"
#include "bng/bng.h"

#include "viz_output_event.h"
#include "scheduler.h"
#include "world.h"

using namespace std;
using namespace DMUtil;
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


void BngDataToDatamodelConverter::to_data_model(const World* world_, Value& mcell_node) {

  reset();

  world = world_;
  bng_engine = &(world->bng_engine);

  // similarly as in pymcell converter, maximum effort is given to conversion and
  // if some parts won't go through, we are trying to generate
  CHECK(convert_species(mcell_node));
  CHECK(convert_rxns(mcell_node));

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


void BngDataToDatamodelConverter::convert_species(Value& mcell_node) {
  Value& define_molecules = mcell_node[KEY_DEFINE_MOLECULES];
  json_add_version(define_molecules, JSON_DM_VERSION_1638);
  Value& molecule_list = define_molecules[KEY_MOLECULE_LIST];

  for (const BNG::Species& s: bng_engine->get_all_species().get_species_vector()) {
    // ALL_* species re ignored and reactive surfaces are converted as surface classes
    if (is_species_superclass(s.name) || s.is_reactive_surface()) {
      continue;
    }
    Value species_node;
    CHECK(convert_single_species(s, species_node));
    molecule_list.append(species_node);
  }
}


void BngDataToDatamodelConverter::convert_single_species(const BNG::Species& s, Value& species_node) {

  json_add_version(species_node, JSON_DM_VERSION_1632);

  Value& display = species_node[KEY_DISPLAY];
  display[KEY_EMIT] = 0.0;
  Value& color_value = display[KEY_COLOR];
  Vec3 c = get_next_color();
  color_value.append(c.r);
  color_value.append(c.g);
  color_value.append(c.b);
  display[KEY_GLYPH] = VALUE_GLYPH_SPHERE_1;
  display[KEY_SCALE] = 1;

  species_node[KEY_BNGL_COMPONENT_LIST] = Value(Json::arrayValue); // empty array;
  species_node[KEY_MOL_BNGL_LABEL] = "";
  species_node[KEY_DESCRIPTION] = "";
  species_node[KEY_MOL_TYPE] = s.is_vol() ? VALUE_MOL_TYPE_3D : VALUE_MOL_TYPE_2D;

  const VizOutputEvent* viz =
      dynamic_cast<const VizOutputEvent*>(world->scheduler.find_next_event_with_type_index(EVENT_TYPE_INDEX_VIZ_OUTPUT));
  if (viz == nullptr || viz->should_visualize_all_species()) {
    // no visualization/or enabled globally
    species_node[KEY_EXPORT_VIZ] = false;
  }
  else {
    species_node[KEY_EXPORT_VIZ] = viz->species_ids_to_visualize.count(s.species_id) == 1;
  }

  species_node[KEY_CUSTOM_SPACE_STEP] = "";
  species_node[KEY_MAXIMUM_STEP_LENGTH] = "";
  species_node[KEY_TARGET_ONLY] = false;
  species_node[KEY_DIFFUSION_CONSTANT] = DMUtil::f_to_string(s.D);
  species_node[KEY_SPATIAL_STRUCTURE] = "None";
  species_node[KEY_MOL_NAME] = s.name;
  species_node[KEY_CUSTOM_TIME_STEP] = "";
}


void BngDataToDatamodelConverter::convert_single_rxn_rule(const BNG::RxnRule& r, Value& rxn_node) {
  assert(r.type == BNG::RxnType::Standard);
  json_add_version(rxn_node, JSON_DM_VERSION_1330);

  // name is put into rnx_name and name is the string of the reaction
  rxn_node[KEY_RXN_NAME] = (r.name != NAME_NOT_SET) ? r.name : "";

  // LATER: maybe find an opposite reaction and generate it as reversible
  rxn_node[KEY_RXN_TYPE] = VALUE_IRREVERSIBLE;

  string reactants = r.reactants_to_str(bng_engine->get_data());
  rxn_node[KEY_REACTANTS] = reactants;
  string products = r.products_to_str(bng_engine->get_data());
  if (products == "") {
    products = VALUE_NULL;
  }
  rxn_node[KEY_PRODUCTS] = products;
  rxn_node[KEY_NAME] = reactants + " -> " + products;

  if (r.variable_rates.empty()) {
    rxn_node[KEY_FWD_RATE] = DMUtil::f_to_string(r.rate_constant);
    rxn_node[KEY_BKWD_RATE] = "";

    CHECK_PROPERTY_W_NAME(r.variable_rates.empty(), r.name);
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

    rxn_node[KEY_VARIABLE_RATE] = "var_rate_" + DMUtil::to_id(r.name) + "_" + to_string(rxn_counter);
    rxn_counter++;

    stringstream text;
    // initial value for time 0
    text << 0.0 << "\t" << r.rate_constant << "\n";
    for (const BNG::RxnRateInfo& ri: r.variable_rates) {
      text << ri.time * world->config.time_unit  << "\t" << ri.rate_constant << "\n";
    }
    rxn_node[KEY_VARIABLE_RATE_TEXT] = text.str();
  }
}


string BngDataToDatamodelConverter::get_surface_class_name(const BNG::RxnRule& r) {
  const string& reactant_name = bng_engine->get_data().get_molecule_type(r.reactants[0].get_simple_species_mol_type_id()).name;
  return DMUtil::remove_surf_class_prefix(reactant_name, r.name);
}


void BngDataToDatamodelConverter::convert_single_surface_class(const BNG::RxnRule& base_rxn, Value& surface_class) {
  surface_class[KEY_DESCRIPTION] = "";
  json_add_version(surface_class, JSON_DM_VERSION_1330);

  string name = get_surface_class_name(base_rxn);
  assert(processed_surface_classes.count(name) == 0 && "This surface class was already processed");
  processed_surface_classes.insert(name);
  surface_class[KEY_NAME] = name;

  Value& surface_class_prop_list = surface_class[KEY_SURFACE_CLASS_PROP_LIST];

  // find all rxn rules with the same name
  vector<const BNG::RxnRule*> rxns_belonging_to_surf_class;
  for (const BNG::RxnRule* rxn_rule: bng_engine->get_all_rxns().get_rxn_rules_vector()) {
    if (get_surface_class_name(*rxn_rule) == name) {
      rxns_belonging_to_surf_class.push_back(rxn_rule);
    }
  }

  // and convert them all as properties of this surface class
  for (const BNG::RxnRule* rxn_rule: rxns_belonging_to_surf_class) {
    Value sc_item;
    json_add_version(sc_item, JSON_DM_VERSION_1756);
    CONVERSION_CHECK(rxn_rule->reactants[0].is_simple(), "Surface class reactant must be simple for now.");

    const string& reactant_name = bng_engine->get_data().get_molecule_type(rxn_rule->reactants[0].get_simple_species_mol_type_id()).name;
    if (is_species_superclass(reactant_name)) {
      sc_item[KEY_AFFECTED_MOLS] = reactant_name;
      sc_item[KEY_MOLECULE] = "";
    }
    else {
      sc_item[KEY_AFFECTED_MOLS] = VALUE_SINGLE;
      sc_item[KEY_MOLECULE] = reactant_name;
    }

    sc_item[KEY_SURF_CLASS_ORIENT] = DMUtil::orientation_to_str(rxn_rule->reactants[0].get_orientation());

    switch (rxn_rule->type) {
      case BNG::RxnType::Transparent:
        sc_item[KEY_SURF_CLASS_TYPE] = VALUE_TRANSPARENT;
        break;
      case BNG::RxnType::Reflect:
        sc_item[KEY_SURF_CLASS_TYPE] = VALUE_REFLECTIVE;
        break;
      case BNG::RxnType::Standard:
        if (rxn_rule->is_absorptive_region_rxn()) {
          sc_item[KEY_SURF_CLASS_TYPE] = VALUE_ABSORPTIVE;
        }
        else {
          assert(false);
        }
        break;
      default:
        CONVERSION_UNSUPPORTED("Unexpected reaction type");
    }

    sc_item[KEY_CLAMP_VALUE] = "0";
    sc_item[KEY_NAME] = ""; // name is ignored by the datamodel to mdl converter anyway

    surface_class_prop_list.append(sc_item);
  }
}


void BngDataToDatamodelConverter::convert_rxns(Value& mcell_node) {

  // --- define_reactions --- (this section is always present)
  Value& define_reactions = mcell_node[KEY_DEFINE_REACTIONS];
  json_add_version(define_reactions, JSON_DM_VERSION_1638);
  Value& reaction_list = define_reactions[KEY_REACTION_LIST];
  reaction_list = Value(Json::arrayValue);

  // --- define_surface_classes --- (only when needed)
  Value& define_surface_classes = mcell_node[KEY_DEFINE_SURFACE_CLASSES];
  json_add_version(define_surface_classes, JSON_DM_VERSION_1638);
  Value& surface_class_list = define_surface_classes[KEY_SURFACE_CLASS_LIST];
  surface_class_list = Value(Json::arrayValue);

  for (const BNG::RxnRule* rxn_rule: bng_engine->get_all_rxns().get_rxn_rules_vector()) {
    if (rxn_rule->type == BNG::RxnType::Standard) {
      // this might be a special case of absorptive reaction
      if (!rxn_rule->is_absorptive_region_rxn()) {
        Value rxn_node;
        CHECK(convert_single_rxn_rule(*rxn_rule, rxn_node));
        reaction_list.append(rxn_node);
      }
      else {
        if (processed_surface_classes.count(get_surface_class_name(*rxn_rule)) == 0) {
          Value surface_class;
          CHECK(convert_single_surface_class(*rxn_rule, surface_class));
          surface_class_list.append(surface_class);
        }
      }
    }
    else if (rxn_rule->type == BNG::RxnType::Transparent || rxn_rule->type == BNG::RxnType::Reflect) {
      // this rxn rule defines a surface class
      if (processed_surface_classes.count(get_surface_class_name(*rxn_rule)) == 0) {
        Value surface_class;
        CHECK(convert_single_surface_class(*rxn_rule, surface_class));
        surface_class_list.append(surface_class);
      }
    }
    else {
      CONVERSION_UNSUPPORTED("AbsorbRegionBorder surf classes are not supported yet");
    }
  }
}


} // namespace MCell

