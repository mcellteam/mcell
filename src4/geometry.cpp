/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <iostream>

#include "bng/bng.h"

#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs_shared.h"
#include "logging.h"
#include "wall_util.h"

#include "partition.h"
#include "geometry.h"
#include "release_event.h"
#include "datamodel_defines.h"

#include "geometry_utils.h"
#include "geometry_utils.inl" // uses get_wall_bounding_box, maybe not include this file
#include "collision_utils.inl"
#include "dump_state.h"

using namespace std;

namespace MCell {

/***************************************************************************
compatible_edges:
  In: array of pointers to walls
      index of first wall
      index of edge in first wall
      index of second wall
      index of edge in second wall
      indices may be -1
  Out: 1 if the edge joins the two walls
       0 if not (i.e. the wall doesn't contain the edge or the edge is
       traversed in the same direction in each or the two walls are
       actually the same wall)
***************************************************************************/
static int compatible_edges(
    const Partition& p,
    const vector<wall_index_t>& faces, int wA, int eA, int wB, int eB) {

  const Vec3 *vA0, *vA1, *vA2, *vB0, *vB1, *vB2;

  const Wall& wall_a = p.get_wall(faces[wA]);
  const Wall& wall_b = p.get_wall(faces[wB]);

  if ((wA < 0) || (eA < 0) || (wB < 0) || (eB < 0))
    return 0;

  vA0 = &p.get_wall_vertex(wall_a, eA);
  if (eA == 2) {
    vA1 = &p.get_wall_vertex(wall_a, 0);
  }
  else {
    vA1 = &p.get_wall_vertex(wall_a, eA + 1);
  }

  if (eA == 0) {
    vA2 = &p.get_wall_vertex(wall_a, 2);
  }
  else {
    vA2 = &p.get_wall_vertex(wall_a, eA - 1);
  }

  vB0 = &p.get_wall_vertex(wall_b, eB);
  if (eB == 2) {
    vB1 = &p.get_wall_vertex(wall_b, 0);
  }
  else {
    vB1 = &p.get_wall_vertex(wall_b, eB + 1);
  }

  if (eB == 0) {
    vB2 = &p.get_wall_vertex(wall_b, 2);
  }
  else {
    vB2 = &p.get_wall_vertex(wall_b, eB - 1);
  }

  return ((vA0 == vB1 && vA1 == vB0 && vA2 != vB2) ||
          (vA0->x == vB1->x && vA0->y == vB1->y && vA0->z == vB1->z &&
           vA1->x == vB0->x && vA1->y == vB0->y && vA1->z == vB0->z &&
           !(vA2->x == vB2->x && vA2->y == vB2->y && vA2->z == vB2->z)));
}


/*****************************************************************************
 have_common_region checks if wall1 and wall2 located on the (same) object
 are part of a common region or not
******************************************************************************/
static bool have_common_region(
    const Partition& p, const GeometryObject& obj, wall_index_t wall1, wall_index_t wall2) {

  const Wall& w1 = p.get_wall(obj.wall_indices[wall1]);
  const Wall& w2 = p.get_wall(obj.wall_indices[wall2]);

  for (region_index_t index: w1.regions) {
    if (w2.regions.count(index) == 1) {
      return true;
    }
  }
  return false;
}


/***************************************************************************
refine_edge_pairs:
  In: the head of a linked list of shared edges
      array of pointers to walls
  Out: No return value.  The best-matching pair of edges percolates up
       to be first in the list.  "Best-matching" means that the edge
       is traversed in different directions by each face, and that the
       normals of the two faces are as divergent as possible.
***************************************************************************/
static void refine_edge_pairs(
    const Partition& p, const GeometryObject& obj, poly_edge *pe, const vector<wall_index_t>& faces) {

#define TSWAP(x, y) temp = (x); (x) = (y); (y) = temp

  int temp;

  pos_t best_align = 2;
  bool share_region = false;
  poly_edge* best_p1 = pe;
  poly_edge* best_p2 = pe;
  int best_n1 = 1;
  int best_n2 = 2;

  poly_edge* p1 = pe;
  int n1 = 1;
  while (p1 != NULL && p1->n >= n1) {
    int wA, eA;
    if (n1 == 1) {
      wA = p1->face[0];
      eA = p1->edge[0];
    }
    else {
      wA = p1->face[1];
      eA = p1->edge[1];
    }

    poly_edge* p2;
    int n2;
    if (n1 == 1) {
      n2 = n1 + 1;
      p2 = p1;
    }
    else {
      n2 = 1;
      p2 = p1->next;
    }

    while (p2 != NULL && p2->n >= n2) {
      int wB, eB;
      if (n2 == 1) {
        wB = p2->face[0];
        eB = p2->edge[0];
      } else {
        wB = p2->face[1];
        eB = p2->edge[1];
      }

      // as soon as we hit an incompatible edge we can break out of the p2 loop
      // and continue scanning the next p1
      if (compatible_edges(p, faces, wA, eA, wB, eB)) {
        const Wall& wall_a = p.get_wall(wA);
        const Wall& wall_b = p.get_wall(wB);
        assert(wall_a.wall_constants_initialized);
        assert(wall_b.wall_constants_initialized);

        pos_t align = dot(wall_a.normal, wall_b.normal);

        // as soon as two walls have a common region we only consider walls who
        // share (any) region. We need to reset the best_align to make sure we
        // don't pick any wall that don't share a region discovered previously
        bool common_region = have_common_region(p, obj, wA, wB);
        if (common_region) {
          if (!share_region) {
            best_align = 2;
          }
          share_region = true;
        }

        if (common_region || !share_region) {
          if (align < best_align) {
            best_p1 = p1;
            best_p2 = p2;
            best_n1 = n1;
            best_n2 = n2;
            best_align = align;
          }
        }
      }
      else {
        break;
      }

      if (n2 == 1) {
        n2++;
      }
      else {
        p2 = p2->next;
        n2 = 1;
      }
    }

    if (n1 == 1) {
      n1++;
    }
    else {
      p1 = p1->next;
      n1 = 1;
    }
  }

  /* swap best match into top spot */
  if (best_align > 1) {
    return; /* No good pairs. */
  }

  TSWAP(best_p1->face[best_n1-1], pe->face[0]);
  TSWAP(best_p1->edge[best_n1-1], pe->edge[0]);
  TSWAP(best_p2->face[best_n2-1], pe->face[1]);
  TSWAP(best_p2->edge[best_n2-1], pe->edge[1]);

#undef TSWAP
}


/***************************************************************************
surface_net:
  In: array of pointers to walls
      integer length of array
  Out: -1 if the surface is a manifold, 0 if it is not, 1 on malloc failure
       Walls end up connected across their edges.
  Note: Two edges must have their vertices listed in opposite order (i.e.
        connect two faces pointing the same way) to be linked.  If more than
        two faces share the same edge and can be linked, the faces with
        normals closest to each other will be linked.  We do not assume that
        the object is connected.  All pieces must be a manifold, however,
        for the entire object to be a manifold.  (That is, there must not
        be any free edges anywhere.)  It is possible to build weird, twisty
        self-intersecting things.  The behavior of these things during a
        simulation is not guaranteed to be well-defined.
***************************************************************************/
static int surface_net(Partition& p, GeometryObject& obj) {

  uint nfaces = obj.wall_indices.size();
  vector<wall_index_t> facelist = obj.wall_indices;

  Edge *e;
  int is_closed = 1;

  struct edge_hashtable eht;
  int nkeys = (3 * nfaces) / 2;
  if (ehtable_init(&eht, nkeys))
    return 1;

  for (uint face_index = 0; face_index < nfaces; face_index++) {

    Wall& w = p.get_wall(facelist[face_index]);

    int k;
    int nedge = 3;

    for (int j = 0; j < nedge; j++) {

      if (j + 1 < nedge)
        k = j + 1;
      else
        k = 0;

      poly_edge pe;
      const Vec3& vert_j = p.get_wall_vertex(w, j);
      pe.v1x = vert_j.x;
      pe.v1y = vert_j.y;
      pe.v1z = vert_j.z;
      const Vec3& vert_k = p.get_wall_vertex(w, k);
      pe.v2x = vert_k.x;
      pe.v2y = vert_k.y;
      pe.v2z = vert_k.z;
      pe.face[0] = face_index;
      pe.edge[0] = j;

      if (ehtable_add(&eht, &pe))
        return 1;
    }
  }

  for (int i = 0; i < nkeys; i++) {
    poly_edge *pep = (eht.data + i);
#ifdef DEBUG_GEOM_OBJ_INITIALIZATION
    dump_poly_edge(i, pep);
#endif
    while (pep != NULL) {
      if (pep->n > 2) {
        refine_edge_pairs(p, obj, pep, facelist);
      }
      if (pep->n >= 2) {
        if (pep->face[0] != -1 && pep->face[1] != -1) {
          if (compatible_edges(p, facelist, pep->face[0], pep->edge[0], pep->face[1], pep->edge[1])) {

            Wall& face0 = p.get_wall(facelist[pep->face[0]]);
            Wall& face1 = p.get_wall(facelist[pep->face[1]]);

            assert(pep->face[0] < (int)facelist.size());
            assert(pep->face[1] < (int)facelist.size());
            face0.nb_walls[pep->edge[0]] = facelist[pep->face[1]];
            face1.nb_walls[pep->edge[1]] = facelist[pep->face[0]];

            Edge e;
            e.forward_index = facelist[pep->face[0]];
            e.backward_index = facelist[pep->face[1]];
            e.edge_num_used_for_init = pep->edge[0];

            assert(face0.wall_constants_initialized);
            assert(face1.wall_constants_initialized);
            e.reinit_edge_constants(p);

            face0.edges[pep->edge[0]] = e;
            face1.edges[pep->edge[1]] = e;
          }

        }
        else {
          is_closed = 0;
        }
      }
      else if (pep->n == 1) {
        is_closed = 0;
        Edge e;

        e.forward_index = facelist[pep->face[0]];
        e.backward_index = WALL_INDEX_INVALID;
        /* Don't call reinit_edge_constants unless both edges are set */
      }
      pep = pep->next;
    }
  }

  ehtable_kill(&eht);
  return -is_closed; /* We use 1 to indicate malloc failure so return 0/-1 */
}


void GeometryObject::initialize_neighboring_walls_and_their_edges(Partition& p) {
  surface_net(p, *this);
}


// returns indices of all vertices in this object that are connected through an edge
const std::set<vertex_index_t>& GeometryObject::get_connected_vertices(
    const Partition& p, const vertex_index_t vi) {

  // check cache first
  auto it = connected_vertices_cache.find(vi);
  if (it != connected_vertices_cache.end()) {
    return it->second;
  }

  // create a new entry
  set<vertex_index_t>& connected_vertices = connected_vertices_cache[vi];

  for (wall_index_t wi: wall_indices) {
    const Wall& w = p.get_wall(wi);

    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      if (w.vertex_indices[i] == vi) {
        // add all other vertices
        for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
          if (w.vertex_indices[k] != vi) {
            connected_vertices.insert(w.vertex_indices[k]);
          }
        }
        break;
      }
    }
  }

  return connected_vertices;
}


wall_index_t GeometryObject::get_wall_for_vertex_pair(
    const Partition& p, const vertex_index_t vi1, const vertex_index_t vi2) {

  assert(vi1 != vi2);

  // check cache first
  UnorderedPair vp(vi1, vi2);
  auto it = vertex_pair_to_wall_cache.find(vp);
  if (it != vertex_pair_to_wall_cache.end()) {
    return it->second;
  }

  // create a new entry
  for (wall_index_t wi: wall_indices) {
    const Wall& w = p.get_wall(wi);

    uint cnt = 0;
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      if (w.vertex_indices[i] == vi1 || w.vertex_indices[i] == vi2) {
        cnt++;
      }
    }
    if (cnt == 2) {
      vertex_pair_to_wall_cache[vp] = w.index;
      return w.index;
    }
  }

  return WALL_INDEX_INVALID;
}


// checks all walls and their regions and if all are only transparent to all molecules,
// sets member is_fully_transparent to true
void GeometryObject::initialize_is_fully_transparent(Partition& p) {

  is_fully_transparent = false;

  const BNG::SpeciesRxnClassesMap* all_molecules_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(p.get_all_species().get_all_molecules_species_id());

  if (all_molecules_rxns == nullptr) {
    // no transparent surface classes for all molecules were defined at all
    return;
  }

  set<species_id_t> transparent_surf_classes_cache;

  // for each wall
  for (wall_index_t wi: wall_indices) {
    const Wall& w = p.get_wall(wi);

    bool is_transparent = false;

    // check all regions
    for (region_index_t ri: w.regions) {
      const Region& reg = p.get_region(ri);
      if (!reg.is_reactive()) {
        continue;
      }

      if (transparent_surf_classes_cache.count(reg.species_id) != 0) {
        is_transparent = true;
        continue;
      }

      auto reactions_it = all_molecules_rxns->find(reg.species_id);
      if (reactions_it == all_molecules_rxns->end()) {
        // no reactions for this type of region -> wall is not transparent
        return;
      }

      BNG::RxnClass* rxn_class = reactions_it->second;

      if (rxn_class->is_transparent_type()) {
        transparent_surf_classes_cache.insert(reg.species_id);
        is_transparent = true;
      }
    }

    if (!is_transparent) {
      // wall has other surface classes -> object is not fully transparent
      return;
    }
  }

  // ok, all walls are fully transparent
  is_fully_transparent = true;
}


void GeometryObject::dump(const Partition& p, const std::string ind) const {
  cout << ind <<
      "GeometryObject: id:" << id << ", name:" << name <<
      ", compartment_id " << vol_compartment_id <<
      ", is_used_in_mol_rxn_counts " << is_used_in_mol_rxn_counts << "\n";

  for (wall_index_t i: wall_indices) {
    cout << ind << "  " << i << ": \n";
    p.get_wall(i).dump(p, ind + "    ");
  }
}


void GeometryObject::dump_array(const Partition& p, const std::vector<GeometryObject>& vec) {
  cout << "GeometryObject array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump(p, "  ");
  }
}


void GeometryObject::to_data_model_as_geometrical_object(
    const Partition& p, const SimulationConfig& config,
    Json::Value& object,
    std::set<rgba_t>& used_colors) const {

  object[KEY_NAME] = DMUtils::remove_obj_name_prefix(parent_name, name);

  bool first = true; // to indicate when to use a comma

  // first note which vertices this object uses
  uint_set<vertex_index_t> used_vertex_indices;
  for (wall_index_t wall_index: wall_indices) {
    const Wall& w = p.get_wall(wall_index);

    for (vertex_index_t vertex_index: w.vertex_indices) {
      used_vertex_indices.insert(vertex_index);
    }
  }

  // then generate vertices and remember mapping
  Json::Value& vertex_list = object[KEY_VERTEX_LIST];
  map<vertex_index_t, uint> map_vertex_index_to_vertex_list_index;
  uint current_index_in_vertex_list = 0;
  for (vertex_index_t i = 0; i < p.get_geometry_vertex_count(); i++) {

    if (used_vertex_indices.count(i) == 1) {
      // define mapping vertex_index -> index in vertex array
      map_vertex_index_to_vertex_list_index[i] = current_index_in_vertex_list;
      current_index_in_vertex_list++;

      // append triple x, y, z
      Vec3 pos = p.get_geometry_vertex(i) * Vec3(p.config.length_unit);
      DMUtils::append_triplet(vertex_list, pos.x, pos.y, pos.z);
    }
  }

  // element connections - they correspond to the ordering of walls
  Json::Value& element_connections = object[KEY_ELEMENT_CONNECTIONS];
  for (wall_index_t wall_index: wall_indices) {
    const Wall& w = p.get_wall(wall_index);

    Json::Value vertex_indices;
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      assert(map_vertex_index_to_vertex_list_index.count(w.vertex_indices[i]) == 1);
      vertex_indices.append(map_vertex_index_to_vertex_list_index[w.vertex_indices[i]]);
    }

    element_connections.append(vertex_indices);
  }

  // surface regions
  insertion_ordered_set<region_index_t> region_indices;
  map<wall_index_t, uint> map_wall_index_to_order_index_in_object;
  for (size_t i = 0; i < wall_indices.size(); i++) {

    map_wall_index_to_order_index_in_object[wall_indices[i]] = i;

    const Wall& w = p.get_wall(wall_indices[i]);
    for (region_index_t region_index: w.regions) {
      region_indices.insert_ordered(region_index);
    }
  }

  vector<Json::Value> surface_regions_vec;
  for (region_index_t region_index: region_indices.get_as_vector()) {
    const Region& reg = p.get_region(region_index);

    if (reg.name_has_suffix_ALL()) {
      continue;
    }

    Json::Value surface_region;

    string name = DMUtils::get_surface_region_name(reg.name);
    surface_region[KEY_NAME] = name;

    Json::Value include_elements;

    for (auto it: reg.walls_and_edges) {
       include_elements.append(map_wall_index_to_order_index_in_object[it.first]);
    }
    surface_region[KEY_INCLUDE_ELEMENTS] = include_elements;
    surface_regions_vec.push_back(surface_region);
  }

  // we may create the define_surface_regions only when there is at least one region
  if (!surface_regions_vec.empty()) {
    Json::Value& define_surface_regions = object[KEY_DEFINE_SURFACE_REGIONS];

    for (auto& surface_region: surface_regions_vec) {
      define_surface_regions.append(surface_region);
    }
  }

  // colors
  if (default_color != DEFAULT_COLOR || !wall_specific_colors.empty()) {
    used_colors.insert(default_color);
    std::map<rgba_t, int> colors_this_object;
    colors_this_object[default_color] = 0;

    // collect all wall_specific_colors and define their ID
    for (const auto& pair_index_color: wall_specific_colors) {
      rgba_t wall_color = pair_index_color.second;
      auto it = colors_this_object.find(wall_color);
      if (it == colors_this_object.end()) {
        // new item
        int local_index = colors_this_object.size();
        colors_this_object[wall_color] = local_index;
      }
    }

    // generate material names for this object, must be sorted by their id
    vector<rgba_t> sorted_colors;
    sorted_colors.resize(colors_this_object.size());
    for (const auto& pair_color_local_index: colors_this_object) {
      sorted_colors[pair_color_local_index.second] = pair_color_local_index.first;
    }
    Json::Value& material_names = object[KEY_MATERIAL_NAMES];
    for (const rgba_t color: sorted_colors) {
      material_names.append(DMUtils::color_to_mat_name(color));
    }

    // and finally assign element_material_indices
    Json::Value& element_material_indices = object[KEY_ELEMENT_MATERIAL_INDICES];
    for (wall_index_t wi: wall_indices) {
      rgba_t color;
      auto it = wall_specific_colors.find(wi);
      if (it != wall_specific_colors.end()) {
        color = it->second;
        used_colors.insert(color);
      }
      else {
        color = default_color;
      }
      assert(colors_this_object.count(color) == 1);
      element_material_indices.append(colors_this_object[color]);
    }
  }
  else {
    // material_names must not be empty
    object[KEY_MATERIAL_NAMES].append(Json::Value(KEY_VALUE_MEMBRANE));
  }
}


void GeometryObject::to_data_model_as_model_object(
    const Partition& p, Json::Value& model_object) const {
  model_object[KEY_DESCRIPTION] = "";
  model_object[KEY_OBJECT_SOURCE] = VALUE_BLENDER;
  model_object[KEY_DYNAMIC_DISPLAY_SOURCE] = "script";
  model_object[KEY_SCRIPT_NAME] = "";

  string obj_name = DMUtils::remove_obj_name_prefix(parent_name, name);
  model_object[KEY_NAME] = obj_name;

  // set defaults that may be overwritten
  model_object[KEY_MEMBRANE_NAME] = "";
  model_object[KEY_PARENT_OBJECT] = "";

  const BNG::BNGData& bng_data = p.get_all_species().get_bng_data();

  if (represents_compartment()) {
    model_object[KEY_IS_BNGL_COMPARTMENT] = true;

    const BNG::Compartment& comp3d = bng_data.get_compartment(vol_compartment_id);
    assert(comp3d.is_3d);
    if (comp3d.name != obj_name) {
      mcell_error(
          "For data model export, name of object %s must be the same as its compartment name %s.",
          obj_name.c_str(), comp3d.name.c_str());
    }

    if (comp3d.has_parent()) {
      const BNG::Compartment& comp_parent = bng_data.get_compartment(comp3d.parent_compartment_id);
      if (!comp_parent.is_3d) {
        // parent of this 3d object is its membrane
        model_object[KEY_MEMBRANE_NAME] = comp_parent.name;

        if (comp_parent.has_parent()) {
          const BNG::Compartment& comp3d_parent =
              bng_data.get_compartment(comp_parent.parent_compartment_id);
          assert(comp3d_parent.is_3d);
          // parent of this 3d object is its membrane
          model_object[KEY_PARENT_OBJECT] = comp3d_parent.name;
        }
      }
      else {
        // parent is this 3d object, membrane has no name
        model_object[KEY_PARENT_OBJECT] = comp_parent.name;
      }
    }
  }
  else {
    model_object[KEY_IS_BNGL_COMPARTMENT] = false;
  }

  model_object[KEY_DYNAMIC] = false;
}


// checks only in debug mode whether the wall index belongs to this object
rgba_t GeometryObject::get_wall_color(const wall_index_t wi) const {
  assert(!wall_indices.empty());
  assert(wi >= wall_indices.front() && wi <= wall_indices.back());

  auto it = wall_specific_colors.find(wi);
  if (it == wall_specific_colors.end()) {
    return default_color;
  }
  else {
    return it->second;
  }
}


// checks only in debug mode whether the wall index belongs to this object
void GeometryObject::set_wall_color(const wall_index_t wi, const rgba_t color) {
  assert(!wall_indices.empty());
  assert(wi >= wall_indices.front() && wi <= wall_indices.back());

  wall_specific_colors[wi] = color;
}


void InitialSurfaceReleases::to_data_model(
    const BNG::SpeciesContainer& all_species,
    Json::Value& initial_region_molecules
) const {
  initial_region_molecules[KEY_MOLECULE] = all_species.get(species_id).name;
  initial_region_molecules[KEY_ORIENT] = DMUtils::orientation_to_str(orientation);
  if (const_num_not_density) {
    initial_region_molecules[KEY_MOLECULE_NUMBER] = to_string(release_num);
  }
  else {
    initial_region_molecules[KEY_MOLECULE_DENSITY] = f_to_str(release_density);
  }
}


void InitialSurfaceReleases::dump(const std::string ind) const {
  cout << ind <<
      "species_id: " << species_id <<
      ", orientation: " << orientation <<
      ", const_num_not_density: " << const_num_not_density <<
      ", release_num: " << release_num <<
      ", release_density: " << release_density <<
      "\n";
}


void Region::init_from_whole_geom_object(const GeometryObject& obj) {
  name = obj.name + REGION_ALL_SUFFIX_W_COMMA;
  geometry_object_id = obj.id;

  // simply add all walls
  for (const wall_index_t& wi: obj.wall_indices) {
    add_wall_to_walls_and_edges(wi, false);
  }
}


void Region::init_surface_region_edges(const Partition& p) {
  assert(!walls_and_edges.empty());

  // must be run after edge must be initializaion,
  // however, not all edges need to be initialized
  for (auto& it: walls_and_edges) {
    const Wall& w = p.get_wall(it.first);
    for (edge_index_t ei = 0; ei < EDGES_IN_TRIANGLE; ei++) {
      const Edge& edge = w.edges[ei];
      // uninitialized edges are considered to be bordering
      if (!edge.is_initialized()) {
        it.second.insert(ei);
      }
      else {
        // is one of the forward or backward walls in our region?
        // if not, then the edge is on the border of this region
        if (walls_and_edges.count(edge.forward_index) == 0 ||
            walls_and_edges.count(edge.backward_index) == 0) {
          it.second.insert(ei);
        }
      }
    }
  }
}


/* tetrahedralVol returns the (signed) volume of the tetrahedron spanned by
 * the vertices a, b, c, and d.
 * The formula was taken from "Computational Geometry" (2nd Ed) by J. O'Rourke
 */
static pos_t tetrahedral_volume(
    const Vec3& a, const Vec3& b,
    const Vec3& c, const Vec3& d
) {
  // TODO: maybe use vector operators
  return (pos_t)1 / (pos_t)6 *
      (-1 * (a.z - d.z) * (b.y - d.y) * (c.x - d.x) +
      (a.y - d.y) * (b.z - d.z) * (c.x - d.x) +
      (a.z - d.z) * (b.x - d.x) * (c.y - d.y) -
      (a.x - d.x) * (b.z - d.z) * (c.y - d.y) -
      (a.y - d.y) * (b.x - d.x) * (c.z - d.z) +
      (a.x - d.x) * (b.y - d.y) * (c.z - d.z));
}


/***************************************************************************
  Based on is_manifold:
  In: r: A region. This region must already be painted on walls. The edges must
         have already been added to the object (i.e. sharpened).
      count_regions_flag: This is usually set, unless we are only checking
                          volumes for dynamic geometries.
  Out: 1 if the region is a manifold, 0 otherwise.
  Note: by "manifold" we mean "orientable compact two-dimensional
        manifold without boundaries embedded in R3"
***************************************************************************/
void Region::initialize_volume_info_if_needed(const Partition& p) {
  if (volume_info_initialized) {
    return;
  }

  if (walls_and_edges.empty()) {
    mcell_internal_error("Region '%s' has NULL wall array!", name.c_str());
  }

  // initialize bounding box
  compute_bounding_box(p, bounding_box_llf, bounding_box_urb);

  // use the center of the region bounding box as reference point for
  // computing the volume
  Vec3 d = (bounding_box_llf + bounding_box_urb) * Vec3(0.5);

  perf() << "Computing volume of region " << name << "\n";

  volume = 0;
  pos_t current_volume = 0;
  for (auto it: walls_and_edges) {

    if (!it.second.empty()) {
      // does this wall represent a region border?
      region_is_manifold = false;
      volume_info_initialized = true;
    }

    const Wall& w = p.get_wall(it.first);

    // compute volume of tetrahedron with w as its face
    current_volume += tetrahedral_volume(
        p.get_wall_vertex(w, 0), p.get_wall_vertex(w, 1), p.get_wall_vertex(w, 2), d);
  }

  perf() << "  - volume computed\n";

  volume = current_volume;
  region_is_manifold = true;
  volume_info_initialized = true;
}


void Region::initialize_wall_subpart_mapping_if_needed(const Partition& p) {
  if (walls_per_subpart_initialized) {
    return;
  }

  perf() << "Preparing wall subpart mapping for region " << name << "\n";

  // first initialize region wall subpart mapping
  walls_per_subpart.clear();

  // FIXME: update when a wall of this region is moved
  for (const auto& wall_and_edges_pair: walls_and_edges) {
    const Wall& w = p.get_wall(wall_and_edges_pair.first);

    for (subpart_index_t si: w.present_in_subparts) {
      walls_per_subpart[si].push_back(w.index);
    }
  }

  perf() << "  - wall subpart mapping prepared\n";

  walls_per_subpart_initialized = true;
}

// returns true if the waypoint with current_waypoint_index is in this region
// similar code as in Partition::create_waypoint
bool Region::initialize_region_waypoint(
    Partition& p,
    const IVec3& current_waypoint_index,
    const bool use_previous_waypoint,
    const IVec3& previous_waypoint_index,
    const bool previous_waypoint_present_in_region
) {

  // waypoint is always in the center of a subpartition
  // thanks to the extra margin on each side, we might be out of the partition
  if (p.is_valid_waypoint_index(current_waypoint_index)) {
    Waypoint& waypoint = p.get_waypoint(current_waypoint_index);


    bool previous_waypoint_use_ok = false;
    uint num_crossed;
    if (use_previous_waypoint) {
      const Waypoint& previous_waypoint = p.get_waypoint(previous_waypoint_index);

      bool must_redo_test = false;

      num_crossed = CollisionUtils::get_num_crossed_region_walls(
          p, waypoint.pos, previous_waypoint.pos,
          *this, must_redo_test
      );

      if (!must_redo_test) {
        previous_waypoint_use_ok = true;
      }
    }

    if (previous_waypoint_use_ok) {
      // ok, test passed safely
      if ((num_crossed % 2 == 0 && previous_waypoint_present_in_region) ||
          (num_crossed % 2 == 1 && !previous_waypoint_present_in_region)) {

        // still in the same region
        waypoints_in_this_region.insert(current_waypoint_index);
        return true;
      }
      else {
        return false;
      }
    }

    bool must_redo_test = false;
    bool inside = false;
    do {
      inside = CollisionUtils::is_point_inside_region_no_waypoints(p, waypoint.pos, *this, must_redo_test);
      if (must_redo_test) {
        // updates values referenced by waypoint
        p.move_waypoint_because_positioned_on_wall(current_waypoint_index);
      }
    } while (must_redo_test);

    if (inside) {
      waypoints_in_this_region.insert(current_waypoint_index);
      return true;
    }
  }

  return false;
}

void Region::initialize_region_waypoints_if_needed(Partition& p) {

  if (region_waypoints_initialized) {
    return;
  }

  perf() << "Initializing waypoints for region " << name << "\n";

  assert(walls_per_subpart_initialized);

  // get initial waypoint
  waypoints_in_this_region.clear();

  IVec3 llf_waypoint_index;
  subpart_index_t subpart_index = p.get_subpart_index(bounding_box_llf);
  p.get_subpart_3d_indices_from_index(subpart_index, llf_waypoint_index);

  // then compute how many waypoints in each dimension we need to check
  Vec3 region_dims = bounding_box_urb - bounding_box_llf;
  // num_waypoints need to be incremented by 2 - for each of the side regions
  IVec3 num_waypoints = region_dims / Vec3(p.config.subpart_edge_length) + Vec3(2);

  bool previous_waypoint_present = false;
  bool use_previous_waypoint = false;
  IVec3 previous_waypoint_index;

  for (int x = 0; x < num_waypoints.x; x++) {
    for (int y = 0; y < num_waypoints.y; y++) {
      for (int z = 0; z < num_waypoints.z; z++) {

        IVec3 current_waypoint_index( llf_waypoint_index + IVec3(x, y, z) );

        previous_waypoint_present = initialize_region_waypoint(
            p, current_waypoint_index,
            use_previous_waypoint, previous_waypoint_index, previous_waypoint_present
        );

        use_previous_waypoint = true;
        previous_waypoint_index = current_waypoint_index;
      }

      // starting a new line - start from scratch
      // NOTE: maybe we do not need to restart this every line
      use_previous_waypoint = false;
      previous_waypoint_present = false;
    }
  }

  region_waypoints_initialized = true;

  perf() << "  - region waypoints initialized\n";
}



bool Region::is_point_inside(Partition& p, const Vec3& pos) {

  initialize_volume_info_if_needed(p);
  initialize_wall_subpart_mapping_if_needed(p);
  initialize_region_waypoints_if_needed(p);

  if (!is_manifold()) {
    mcell_error("Cannot check whether a point is a non-manifold volumetric region %s. (possibly during a molecule release)", name.c_str());
  }

  // get a waypoint close to this position
  IVec3 waypoint_index;
  subpart_index_t subpart_index = p.get_subpart_index(pos); // TODO: get indices directly
  p.get_subpart_3d_indices_from_index(subpart_index, waypoint_index);
  const Waypoint& waypoint = p.get_waypoint(waypoint_index);

  map<geometry_object_index_t, uint> num_crossed_walls_per_object;
  bool must_redo_test = false;
  uint num_crossed = CollisionUtils::get_num_crossed_region_walls(
      p, pos, waypoint.pos,
      *this, must_redo_test
  );
  if (!must_redo_test) {
    // ok, using waypoint passed
    bool waypoint_is_inside_this_region = waypoints_in_this_region.count(waypoint_index) != 0;
    if (waypoint_is_inside_this_region) {
      return num_crossed % 2 == 0;
    }
    else {
      return num_crossed % 2 == 1;
    }
  }
  else {
    // let's try to compute the containment without waypoints as a fallback
    bool inside = CollisionUtils::is_point_inside_region_no_waypoints(p, pos, *this, must_redo_test);
    release_assert(!must_redo_test);
    return inside;
  }
}


void Region::dump(const std::string ind, const bool with_geometry) const {
  cout << ind << "Region : " <<
      "name:" << name <<
      ", species_id: " << ((species_id == SPECIES_ID_INVALID) ? string("invalid") : to_string(species_id)) <<
      "\n";

  cout << ind << "  "  << "initial_region_molecules :\n";
  for (const auto& initial_mols: initial_region_molecules) {
    initial_mols.dump(ind + "    ");
  }

  if (with_geometry) {
    for (auto& wall_it: walls_and_edges) {
      cout << ind << "  " << "wall " << wall_it.first << ", region edges: {";
      for (auto& edge: wall_it.second) {
        cout << edge << ", ";
      }
      cout << "}\n";
    }
  }
}


void Region::dump_array(const std::vector<Region>& vec) {
  cout << "Region array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}


void Region::to_data_model(const Partition& p, Json::Value& modify_surface_region) const {

  if (initial_region_molecules.empty()) {
    DMUtils::add_version(modify_surface_region, VER_DM_2018_01_11_1330);
  }
  else {
    DMUtils::add_version(modify_surface_region, VER_DM_2020_07_12_1600);
  }

  modify_surface_region[KEY_DESCRIPTION] = "";
  modify_surface_region[KEY_OBJECT_NAME] =
      DMUtils::remove_obj_name_prefix(p.get_geometry_object(geometry_object_id).name);

  string region_name = DMUtils::get_region_name(name);
  if (region_name == "ALL") {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_ALL;
    modify_surface_region[KEY_REGION_NAME] = "";
  }
  else {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_SEL;
    modify_surface_region[KEY_REGION_NAME] = region_name;
  }

  modify_surface_region[KEY_NAME] = ""; // don't care

  if (species_id != SPECIES_ID_INVALID) {
    modify_surface_region[KEY_SURF_CLASS_NAME] = p.get_species(species_id).name;
  }

  if (!initial_region_molecules.empty()) {
    Json::Value& initial_region_molecules_list = modify_surface_region[KEY_INITIAL_REGION_MOLECULES_LIST];
    for (const InitialSurfaceReleases& info: initial_region_molecules) {
      Json::Value initial_region_molecules_item;
      info.to_data_model(p.get_all_species(), initial_region_molecules_item);
      initial_region_molecules_list.append(initial_region_molecules_item);
    }
  }
}


namespace Geometry {

double compute_geometry_object_area(const Partition& p, const GeometryObject& obj) {
  double res = 0;
  for (wall_index_t wi: obj.wall_indices) {
    res += p.get_wall(wi).area;
  }

  return res;
}


/***************************************************************************
eval_rel_region_bbox:
  In: release expression for a 3D region release
      place to store LLF corner of the bounding box for the release
      place to store URB corner
  Out: true on success, false on failure.  Bounding box is set based on release
       expression (based boolean intersection of bounding boxes for each
       region).  The function reports failure if any region is unclosed.
***************************************************************************/
// TODO: not checking whether object is closed - use some library for that
bool compute_region_expr_bounding_box(
    World* world, const RegionExprNode* expr,
    Vec3& llf, Vec3& urb
) {
  int count_regions_flag = 1;
  Partition& p = world->get_partition(0);

  if (expr->op == RegionExprOperator::LEAF_SURFACE_REGION) {
    Region& reg = p.get_region_by_id(expr->region_id);
    reg.compute_bounding_box(p, llf, urb);
    return true;
  }
  else if (expr->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    GeometryObject& go = p.get_geometry_object(expr->geometry_object_id);
    Region& reg = p.get_region_by_id(go.encompassing_region_id);
    reg.compute_bounding_box(p, llf, urb);
    return true;
  }
  else {
    Vec3 llf_left, urb_left;
    compute_region_expr_bounding_box(world, expr->left, llf_left, urb_left);

    Vec3 llf_right, urb_right;
    compute_region_expr_bounding_box(world, expr->right, llf_right, urb_right);

    llf = llf_left;
    urb = urb_left;

    if (expr->op == RegionExprOperator::UNION) {
      if (llf.x > llf_right.x) {
        llf.x = llf_right.x;
      }
      if (llf.y > llf_right.y) {
        llf.y = llf_right.y;
      }
      if (llf.z > llf_right.z) {
        llf.z = llf_right.z;
      }
      if (urb.x < urb_right.x) {
        urb.x = urb_right.x;
      }
      if (urb.y < urb_right.y) {
        urb.y = urb_right.y;
      }
      if (urb.z < urb_right.z) {
        urb.z = urb_right.z;
      }
    }
    else if (expr->op == RegionExprOperator::DIFFERENCE) {
      // for difference/subtraction the MCell3 implementation returns
      // the llf_left/urb_left
    }
    else if (expr->op == RegionExprOperator::INTERSECT) {
      if (llf.x < llf_right.x) {
        llf.x = llf_right.x;
      }
      if (llf.y < llf_right.y) {
        llf.y = llf_right.y;
      }
      if (llf.z < llf_right.z) {
        llf.z = llf_right.z;
      }
      if (urb.x > urb_right.x) {
        urb.x = urb_right.x;
      }
      if (urb.y > urb_right.y) {
        urb.y = urb_right.y;
      }
      if (urb.z > urb_right.z) {
        urb.z = urb_right.z;
      }
    }
    else {
      assert(false);
      return false;
    }
    return true;
  }
}



// this is the entry point called from Partition class
void update_moved_walls(
    Partition& p,
    const std::vector<VertexMoveInfo*>& scheduled_vertex_moves,
    // we can compute all the information already from scheduled_vertex_moves,
    // but the keys of the map walls_with_their_moves are the walls that we need to update
    const WallsWithTheirMovesMap& walls_with_their_moves
) {

  // move all vertices
  for (const VertexMoveInfo* move_info: scheduled_vertex_moves) {
    // ignore moves where walls are fixed
    if (!move_info->vertex_walls_are_movable) {
      continue;
    }

    Vec3& vertex_ref = p.get_geometry_vertex(move_info->vertex_index);
    vertex_ref = vertex_ref + move_info->displacement;
    if (! p.in_this_partition(vertex_ref) ) {
      mcell_log("Error: Crossing partitions is not supported yet.\n");
      exit(1);
    }
  }

  // update walls
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    Wall& w = p.get_wall(wall_index);

    // first we need to update all wall constants
    w.initialize_wall_constants(p);

    // then update wall_collision_rejection_data that hold a copy
    p.update_wall_collision_rejection_data(w);

    // reinitialize grid
    if (w.grid.is_initialized()) {
      w.grid.initialize(p, w);
    }
  }


  // edges need to be fixed after all wall have been moved
  // otherwise the edge initialization would be using
  // inconsistent data
  // we need to update also edges of neighboring walls
  uint_set<wall_index_t> walls_to_be_updated;
  // NOTE: we might consider sharing edges in the same way as in MCell3
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    Wall& w = p.get_wall(wall_index);

    walls_to_be_updated.insert(wall_index);
    for (uint n = 0; n < EDGES_IN_TRIANGLE; n++) {
      walls_to_be_updated.insert(w.nb_walls[n]);
    }
  }
  for (wall_index_t wall_index: walls_to_be_updated) {
    if (wall_index == WALL_INDEX_INVALID) {
      continue;
    }
    Wall& w = p.get_wall(wall_index);
    w.initialize_edge_constants(p);
  }
}


void rgba_to_components(
    const rgba_t rgba,
    double& red, double& green, double& blue, double& alpha) {

  const double MAX = 255.0;
  red = (((uint)rgba >> 24) & 0xFF) / MAX;
  green = (((uint)rgba >> 16) & 0xFF) / MAX;
  blue = (((uint)rgba >> 8) & 0xFF) / MAX;
  alpha = ((uint)rgba & 0xFF) / MAX;
}

} /* namespace Geometry */

} /* namespace MCell */
