/******************************************************************************
 *
 * Copyright (C) 2019 by
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

/*
Regex to replace struct member definition by dumping code:

^[ \t]+([^ ]+) ([\*]?)([^ ]+);[^ \t]*(.*)

  cout << ind2 << "\3: \\t\\t" << block->\3 << " [\1\2] \\t\\t\4\\n";

 */

#define DUMP_SCHEDULERS
//#define DUMP_WALLS
//#define DUMP_WAYPOINTS
//#define DUMP_SUBVOLUMES
//#define DUMP_RELEASE_REGION_DATA

#include "dump_state.h"

#include <iostream>
#include <string>
#include <cassert>
#include <string.h>

#include "dyngeom_parse_extras.h"
#include "bng_util.h"

#include "debug_config.h"

#define MAX_ARRAY_ITEMS 16
#define MAX_SUBVOLUMES 27
#define DUMP_ARRAY_NEWLINE_COUNT 8

#define IND_ADD2(ind) (std::string(ind) + "  ").c_str()
#define DECL_IND2(ind) std::string inds = ind; inds += "  ";  const char* ind2 = inds.c_str();

using namespace std;

void dump_species(species* spec, const char* name, const char* comment, const char* ind);
void dump_species_list(int n_species, const char* num_name, species** species_list, const char* name, const char* comment, const char* ind);
void dump_species_item(species* spec, const char* ind);
void dump_object_list(geom_object* obj, const char* name, const char* comment, const char* ind);
void dump_wall_list(wall_list* list, const char* ind);
void dump_wall_array(int num, wall** wall_array, const char* ind);
void dump_object(geom_object* o, const char* ind);
void dump_region_list(region_list *rl, const char* regions_name, const char* comment, const char* ind);

const char* str(const char* ptr) {
  if (ptr != nullptr) {
    return ptr;
  }
  else {
    return "0";
  }
}

static const char* get_sym_name(const sym_entry *s) {
  if (s == nullptr) {
    return "NULL SYMBOL";
  }
  else {
    return s->name;
  }
}

std::ostream & operator<<(std::ostream &out, const timeval &a) {
  out << a.tv_sec << "s, " << a.tv_usec << "us";
  return out;
}

std::ostream & operator<<(std::ostream &out, const vector2 &a) {
  out << "(" << a.u << ", " << a.v << ")";
  return out;
}

std::ostream & operator<<(std::ostream &out, const vector3 &a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}

std::ostream & operator<<(std::ostream &out, const periodic_image &a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}

std::ostream & operator<<(std::ostream &out, const int3D &a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}

std::ostream & operator<<(std::ostream &out, const name_orient* a) {
  if (a == nullptr) {
    out << "NULL";
  }
  else {
    for (const name_orient* curr = a; curr != nullptr; curr = curr->next) {
      out << "(" << ((curr->name == nullptr) ? "NULL" : curr->name) << "," << curr->orient << "), ";
    }

  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const viz_mode_t &a) {
  switch (a) {
    case NO_VIZ_MODE: out << "NO_VIZ_MODE"; break;
    case ASCII_MODE: out << "ASCII_MODE"; break;
    case CELLBLENDER_MODE: out << "CELLBLENDER_MODE"; break;
    default: assert(false);
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const output_timer_type_t &a) {
  switch (a) {
    case OUTPUT_BY_STEP: out << "OUTPUT_BY_STEP"; break;
    case OUTPUT_BY_TIME_LIST: out << "OUTPUT_BY_TIME_LIST"; break;
    case OUTPUT_BY_ITERATION_LIST: out << "OUTPUT_BY_ITERATION_LIST"; break;
    default: assert(false);
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const viz_frame_type_t &a) {
  switch (a) {
    case MOL_POS: out << "MOL_POS"; break;
    case MOL_ORIENT: out << "MOL_ORIENT"; break;
    case ALL_MOL_DATA: out << "ALL_MOL_DATA"; break;
    case NUM_FRAME_TYPES: out << "NUM_FRAME_TYPES"; break;
    default: assert(false);
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const object_type_t &a) {
  switch (a) {
    case META_OBJ: out << "META_OBJ"; break;
    case BOX_OBJ: out << "BOX_OBJ"; break;
    case POLY_OBJ: out << "POLY_OBJ"; break;
    case REL_SITE_OBJ: out << "REL_SITE_OBJ"; break;
    case VOXEL_OBJ: out << "VOXEL_OBJ"; break;
    default: assert(false);
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const overwrite_policy_t &a) {
  switch (a) {
    case FILE_UNDEFINED: out << "FILE_UNDEFINED"; break;
    case FILE_OVERWRITE: out << "FILE_OVERWRITE"; break;
    case FILE_SUBSTITUTE: out << "FILE_SUBSTITUTE"; break;
    case FILE_APPEND: out << "FILE_APPEND"; break;
    case FILE_APPEND_HEADER: out << "FILE_APPEND_HEADER"; break;
    case FILE_CREATE: out << "FILE_CREATE"; break;
    default: assert(false);
  }
  return out;
}



std::ostream & operator<<(std::ostream &out, const count_type_t &a) {
  switch (a) {
    case COUNT_UNSET: out << "COUNT_UNSET"; break;
    case COUNT_DBL: out << "COUNT_DBL"; break;
    case COUNT_INT: out << "COUNT_INT"; break;
    case COUNT_TRIG_STRUCT: out << "COUNT_TRIG_STRUCT"; break;
    default: assert(false);
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const num_expr_list* a) {
  out << "(";
  const num_expr_list* p = a;
  while (p != nullptr) {
    out << p->value << ", ";
    p = p->next;
  }
  out << ")";
  return out;
}


std::ostream & operator<<(std::ostream &out, const double a[4][4]) {
  out << "(";
  for (int x = 0; x < 4; x++) {
    out << "(";
    for (int y = 0; y < 4; y++) {
      out << a[x][y];
      if (y != 3)
        out << ",";
    }
    out << ")";
    if (x != 3)
      out << ",";
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const pointer_hash &a) {
    //out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  out << "(pointer_hash - TODO)";
  return out;
}

std::ostream & operator<<(std::ostream &out, const sym_entry *s) {
  if (s != NULL) {
    // TODO: sym_type - see store_sym
    out << "sym_entry(next:" << (void*)s->next << ", sym_type:" << s->sym_type << ", name:" << s->name << ", value:" << (void*)s->value << ", count:" << s->count << ")";
  }
  else {
    out << "sym_entry(NULL)";
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const reaction_flags &a) {

#define DUMP_ITEM(name) #name << "=" << a.name << ", "
  out
    << DUMP_ITEM(vol_vol_reaction_flag)
    << DUMP_ITEM(vol_surf_reaction_flag)
    << DUMP_ITEM(surf_surf_reaction_flag)
    << DUMP_ITEM(vol_wall_reaction_flag)
    << DUMP_ITEM(vol_vol_vol_reaction_flag)
    << DUMP_ITEM(vol_vol_surf_reaction_flag)
    << DUMP_ITEM(vol_surf_surf_reaction_flag)
    << DUMP_ITEM(vol_surf_surf_reaction_flag)
  ;
#undef DUMP_ITEM
  return out;
}

std::ostream & operator<<(std::ostream &out, const output_expression *e);

void dump_oexpr_element(std::ostream &out, const int flags, void* child, bool left) {

  int specific_flags = flags & ((left) ? OEXPR_LEFT_MASK : OEXPR_RIGHT_MASK);

  if (specific_flags == OEXPR_LEFT_INT || specific_flags == OEXPR_RIGHT_INT) {
    out << *(int*)child << "(" << (void*)child << ")";
  }
  else if (specific_flags == OEXPR_LEFT_DBL || specific_flags == OEXPR_RIGHT_DBL) {
    out << *(double*)child;
  }
  else if (specific_flags == OEXPR_LEFT_OEXPR || specific_flags == OEXPR_RIGHT_OEXPR) {
    out << (const output_expression*)child;
  }
  else if (specific_flags == OEXPR_LEFT_CONST || specific_flags == OEXPR_RIGHT_CONST) {
    out << "TODO const";
  }
  else if (specific_flags == OEXPR_LEFT_REQUEST || specific_flags == OEXPR_RIGHT_REQUEST) {
    out << "TODO request";
  }
  else if (specific_flags == OEXPR_LEFT_TRIG || specific_flags == OEXPR_RIGHT_TRIG) {
    out << "TODO trig";
  }
  else {
    out << ""; /*no children*/
  }

  out.flush();
}


std::ostream & operator<<(std::ostream &out, const output_expression *e) {
  // NOTE: the flags for column_head::expr are somehow messed up... everything seems to be an integer
  if (e != NULL) {
    out <<
        "("
        "[" << ( (e->title != nullptr) ? e->title : "") << ", v:" << e->value << "] ";

    dump_oexpr_element(out, e->expr_flags, e->left, true);
    out << e->oper << " ";
    dump_oexpr_element(out, e->expr_flags, e->right, false);
    out << ")";

  }
  else {
    out << "(NULL)";
  }
  return out;
}




#define DUMP_FLAG(f, mask) if (((f) & (mask)) != 0) res += string(#mask) + ", ";

string get_species_flags_string(uint flags) {
  string res;
  DUMP_FLAG(flags, ON_GRID)
  DUMP_FLAG(flags, IS_SURFACE)
  //DUMP_FLAG(flags, NOT_FREE) == ON_GRID | IS_SURFACE
  DUMP_FLAG(flags, TIME_VARY)
  DUMP_FLAG(flags, CAN_VOLVOLVOL)
  DUMP_FLAG(flags, CAN_VOLVOL)
  DUMP_FLAG(flags, CAN_VOLSURF)
  DUMP_FLAG(flags, CAN_VOLWALL)
  DUMP_FLAG(flags, CAN_SURFSURF)
  DUMP_FLAG(flags, CAN_SURFWALL)
  DUMP_FLAG(flags, CAN_VOLVOLSURF)
  DUMP_FLAG(flags, CANT_INITIATE)
  DUMP_FLAG(flags, COUNT_TRIGGER)
  DUMP_FLAG(flags, COUNT_CONTENTS)
  DUMP_FLAG(flags, COUNT_HITS)
  DUMP_FLAG(flags, COUNT_RXNS)
  DUMP_FLAG(flags, COUNT_ENCLOSED)
  DUMP_FLAG(flags, CAN_VOLSURFSURF)
  DUMP_FLAG(flags, CAN_SURFSURFSURF)
  DUMP_FLAG(flags, SET_MAX_STEP_LENGTH)
  DUMP_FLAG(flags, CAN_REGION_BORDER)
  DUMP_FLAG(flags, REGION_PRESENT)
  DUMP_FLAG(flags, EXTERNAL_SPECIES)
  return res;
}

string get_molecule_flags_string(short flags, bool full_dump = true) {
  string res;
  DUMP_FLAG(flags, TYPE_SURF)
  DUMP_FLAG(flags, TYPE_VOL)
  if (full_dump) {
    DUMP_FLAG(flags, ACT_DIFFUSE)
    DUMP_FLAG(flags, ACT_REACT)
    DUMP_FLAG(flags, ACT_NEWBIE)
    DUMP_FLAG(flags, ACT_CHANGE)
  }
  DUMP_FLAG(flags, ACT_CLAMPED)
  if (full_dump) {
    DUMP_FLAG(flags, IN_SCHEDULE)
    DUMP_FLAG(flags, IN_SURFACE)
  }
  DUMP_FLAG(flags, IN_VOLUME)
  return res;
}

string get_report_type_flags_string(char report_type) {

  string res;
#define DUMP_MASKED_VALUE(f, val) if ( ((f) & REPORT_TYPE_MASK) == val) res += string(#val) + ", ";


  DUMP_MASKED_VALUE(report_type, REPORT_NOTHING)
  DUMP_MASKED_VALUE(report_type, REPORT_CONTENTS)
  DUMP_MASKED_VALUE(report_type, REPORT_RXNS)
  DUMP_MASKED_VALUE(report_type, REPORT_FRONT_HITS)
  DUMP_MASKED_VALUE(report_type, REPORT_BACK_HITS)
  DUMP_MASKED_VALUE(report_type, REPORT_FRONT_CROSSINGS)
  DUMP_MASKED_VALUE(report_type, REPORT_BACK_CROSSINGS)
  DUMP_MASKED_VALUE(report_type, REPORT_ALL_HITS)
  DUMP_MASKED_VALUE(report_type, REPORT_ALL_CROSSINGS)
  DUMP_MASKED_VALUE(report_type, REPORT_CONCENTRATION)
  DUMP_MASKED_VALUE(report_type, REPORT_ELAPSED_TIME)

  DUMP_FLAG(report_type, REPORT_WORLD)
  DUMP_FLAG(report_type, REPORT_ENCLOSED)
  DUMP_FLAG(report_type, REPORT_TRIGGER)

  return res;

}


string get_region_expression_flags_string(char op) {
  string res;
  DUMP_FLAG(op, REXP_NO_OP)
  DUMP_FLAG(op, REXP_UNION)
  DUMP_FLAG(op, REXP_INTERSECTION)
  DUMP_FLAG(op, REXP_SUBTRACTION)
  DUMP_FLAG(op, REXP_LEFT_REGION)
  DUMP_FLAG(op, REXP_RIGHT_REGION)

  return res;
}

#undef DUMP_FLAG

string get_release_shape_name(int8_t release_shape) {
  string res;
#define CASE_ITEM(n) case n: res = #n; break;
  switch (release_shape) {
    CASE_ITEM(SHAPE_UNDEFINED)
    CASE_ITEM(SHAPE_SPHERICAL)
    CASE_ITEM(SHAPE_CUBIC)
    CASE_ITEM(SHAPE_ELLIPTIC)
    CASE_ITEM(SHAPE_RECTANGULAR)
    CASE_ITEM(SHAPE_SPHERICAL_SHELL)
    CASE_ITEM(SHAPE_REGION)
    CASE_ITEM(SHAPE_LIST)
    default: res = "unknown!"; break;
  }
#undef CASE_ITEM
  return res;
}


string get_release_number_method_string(release_number_type_t rt) {
  string res;
#define CASE_ITEM(n) case n: res = #n; break;
  switch (rt) {
    CASE_ITEM(CONSTNUM)
    CASE_ITEM(GAUSSNUM)
    CASE_ITEM(VOLNUM)
    CASE_ITEM(CCNNUM)
    CASE_ITEM(DENSITYNUM)
    default: res = "unknown!"; break;
  }
#undef CASE_ITEM
  return res;
}


void dump_double_array(int num, const char* num_name, double* values, const char* values_name, const char* comment, const char* ind, const int max = MAX_ARRAY_ITEMS) {
  cout << ind << values_name << "[" << num_name << "]: \t\t" << values << "[" << num << "]" << " [double[]] \t\t" << comment;
  for (int i = 0; i < num && i < max; i++) {
    if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
      cout << "\n" << ind << "  ";
    }
    cout << i << ":" << values[i] << ", ";
  }
  if (num >= max) {
    cout << "...";
  }
  cout << "\n";
}


void dump_int_array(int num, const char* num_name, int* values, const char* values_name, const char* comment, const char* ind) {
  cout << ind << values_name << "[" << num_name << "]: \t\t" << values << "[" << num << "]" << " [int[]] \t\t" << comment;
  for (int i = 0; i < num && i < MAX_ARRAY_ITEMS; i++) {
    if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
      cout << "\n" << ind << "  ";
    }
    cout << i << ":" << values[i] << ", ";
  }
  if (num >= MAX_ARRAY_ITEMS) {
    cout << "...";
  }
  cout << "\n";
}


void dump_string_array(int num, const char* num_name, const char** values, const char* values_name, const char* comment, const char* ind) {
  cout << ind << values_name << "[" << num_name << "]: \t\t" << values << "[" << num << "]" << " [char*[]] \t\t" << comment;
  for (int i = 0; i < num && i < MAX_ARRAY_ITEMS; i++) {
    if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
      cout << "\n" << ind << "  ";
    }
    cout << i << ":" << str(values[i]) << ", ";
  }
  if (num >= MAX_ARRAY_ITEMS) {
    cout << "...";
  }
  cout << "\n";
}

void dump_vector3_array(int num, const char* num_name, vector3* values, const char* values_name, const char* comment, const char* ind) {
  cout << ind << values_name << "[" << num_name << "]: \t\t" << values << "[" << num << "]" << " [vector3[]] \t\t" << comment;
  for (int i = 0; i < num && i < MAX_ARRAY_ITEMS; i++) {
    if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
      cout << "\n" << ind << "  ";
    }
    cout << i << ":" << values[i] << ", ";
  }
  if (num >= MAX_ARRAY_ITEMS) {
    cout << "...";
  }
  cout << "\n";
}


void dump_object(geom_object* o, const char* ind) {
  DECL_IND2(ind);

  if (o == nullptr) {
    return;
  }

  cout << ind << "next: *\t\t" << o->next << " [object] \t\t/* Next sibling object */\n";
  cout << ind << "parent: *\t\t" << o->parent << " [object] \t\t/* Parent meta object */\n";
  cout << ind << "first_child: *\t\t" << o->first_child << " [object] \t\t/* First child object */\n";
  cout << ind << "last_child: *\t\t" << o->last_child << " [object] \t\t/* Last child object */\n";
  cout << ind << "sym: *\t\t" << o->sym << " [sym_entry] \t\t/* Symbol hash table entry for this object */\n";
  if (o->last_name != nullptr) {
    cout << ind << "last_name: *\t\t" << o->last_name << " [char] \t\t/* Name of object without pre-pended parent object name */\n";
  }
  else {
    cout << ind << "last_name: *\t\t" << (void*)o->last_name << " [char] \t\t/* Name of object without pre-pended parent object name */\n";
  }
  cout << ind << "object_type: \t\t" << o->object_type << " [object_type_flags_t] \t\t/* Object Type Flags */\n";
  cout << ind << "contents: *\t\t" << o->contents << " [void] \t\t/* Actual physical object, cast according to object_type */\n";
  cout << ind << "num_regions: \t\t" << o->num_regions << " [u_int] \t\t/* Number of regions defined on object */\n";

  //cout << ind << "regions: *\t\t" << o->regions << " [region_list] \t\t/* List of regions for this object */\n";
  dump_region_list(o->regions, "regions", "/* List of regions for this object */", ind2);

  cout << ind << "n_walls: \t\t" << o->n_walls << " [int] \t\t/* Total number of walls in object */\n";
  cout << ind << "n_walls_actual: \t\t" << o->n_walls_actual << " [int] \t\t/* Number of non-null walls in object */\n";
  cout << ind << "n_verts: \t\t" << o->n_verts << " [int] \t\t/* Total number of vertices in object */\n";
  cout << ind << "vertices: **\t\t" << o->vertices << " [vector3] \t\t/* Array of pointers to vertices (linked to all_vertices array) */\n";
  cout << ind << "total_area: \t\t" << o->total_area << " [double] \t\t/* Area of object in length units */\n";
  cout << ind << "n_tiles: \t\t" << o->n_tiles << " [u_int] \t\t/* Number of surface grid tiles on object */\n";
  cout << ind << "n_occupied_tiles: \t\t" << o->n_occupied_tiles << " [u_int] \t\t/* Number of occupied tiles on object */\n";
  cout << ind << "t_matrix: *\t\t" << o->t_matrix << " [double[4][4]] \t\t /* Transformation matrix for object */\n";
  cout << ind << "is_closed: \t\t" << o->is_closed << " [short] \t\t/* Flag that describes the geometry of the polygon object (e.g. for sphere is_closed = 1 and for plane is 0) */\n";
  cout << ind << "periodic_x: \t\t" << o->periodic_x << " [bool] \t\t// This flag only applies to box objects BOX_OBJ. If set\n";
  cout << ind << "periodic_y: \t\t" << o->periodic_y << " [bool] \t\t// any volume molecules encountering the box surface in the x,\n";
  cout << ind << "periodic_z: \t\t" << o->periodic_z << " [bool] \t\t// y or z direction are reflected back into the box as if they  had entered the adjacent neighboring box */\n";

  dump_object_list(o->first_child, "first_child", "", ind2);

#ifdef DUMP_WALLS
  cout << ind << "walls: *\t\t" << o->walls << " [wall] \t\t/* Array of walls in object */\n";
  cout << ind << "wall_p: **\t\t" << o->wall_p << " [wall] \t\t// Array of ptrs to walls in object (used at run-time)\n";
  dump_wall_array(o->n_walls, o->wall_p, ind2);
#endif
}


void dump_object_list(geom_object* obj, const char* name, const char* comment, const char* ind) {
  cout << ind << name << " " << comment << "\n";
  geom_object* curr = obj;
  int i = 0;
  while (curr != NULL) {
    cout << ind << i << ":\n";
    dump_object(curr, IND_ADD2(ind));
    i++;
    curr = curr->next;
  }
}


void dump_edge(edge* e, const char* ind, const bool for_diff) {
  DECL_IND2(ind);


  if (e == nullptr) {
    cout << ind << "NULL\n";
  }

  if (for_diff) {
    cout << ind <<
      "Edge: translate: " << e->translate <<
      ", cos_theta: " << e->cos_theta <<
      ", sin_theta: " << e->sin_theta <<
      ", edge_num_used_for_init: " << e->edge_num_used_for_init << "\n";
  }
  else {

    cout << ind2 << "translate: \t\t" << e->translate << " [vector2] \t\t /* Translation vector between coordinate systems */\n";
    cout << ind2 << "cos_theta: \t\t" << e->cos_theta << " [double] \t\t         /* Cosine of angle between coordinate systems */\n";
    cout << ind2 << "sin_theta: \t\t" << e->sin_theta << " [double] \t\t         /* Sine of angle between coordinate systems */\n";

    cout << ind2 << "length: \t\t" << e->length << " [double] \t\t   /* Length of the shared edge */\n";
    cout << ind2 << "length_1: \t\t" << e->length_1 << " [double] \t\t /* Reciprocal of length of shared edge */\n";
    cout << ind2 << "edge_num_used_for_init: \t\t" << e->edge_num_used_for_init << " [int] \t\t // edge used for initialization\n";
  }

}

void dump_wall(wall* w, const char* ind, const bool for_diff) {

  if (for_diff) {
    if (w == nullptr) {
      cout << "wall: NULL!!!";
      return;
    }
    cout << "wall[side: " << w->side << "]: ";
    for (uint i = 0; i < 3; i++) {
      cout << *w->vert[i];
      if (i != 2) {
        cout << ", ";
      }
    }
    cout << "\n";
  }
  else {
    if (w == nullptr) {
      return;
    }

    DECL_IND2(ind);

    cout << ind << "this wall: " << (void*)w << "\n";
    cout << ind << "next: *\t\t" << (void*)w->next << " [wall] \t\t/* Next wall in the universe */\n";

    cout << ind << "surf_class_head: *\t\t" << (void*)w->surf_class_head << " [surf_class_list] \t\t/* linked list of surface classes for this wall (multiple surface classes may come from the overlapping regions */\n";
    cout << ind << "num_surf_classes: \t\t" << w->num_surf_classes << " [int] \t\t/* number of attached surface classes */\n";

    cout << ind << "side: \t\t" << w->side << " [int] \t\t/* index of this wall in its parent object */\n";

    if (w->vert != NULL) {
      cout << ind << "vert[3]: 0: " << *w->vert[0] << ", 1: " << *w->vert[1] << ", 2: " << *w->vert[2] << "\n";
    }
    else {
      cout << ind << "vert: \t\t" << (void*)w->vert << " [vector[3]] \t\t/* Array of pointers to vertices */\n";
    }


    cout << ind << "uv_vert1_u: \t\t" << w->uv_vert1_u << " [double] \t\t/* Surface u-coord of 2nd corner (v=0) */\n";
    cout << ind << "uv_vert2: \t\t" << w->uv_vert2 << " [vector2] \t\t/* Surface coords of third corner */\n";

    cout << ind << "edges[3]: *\t\t" << (void*)w->edges << " [*edges[3]] \t\t/*  /* Array of pointers to each edge. */ // TODO */\n";
    for (int i = 0; i < 3; i++) {
      cout << ind << i << ":\n";
      dump_edge(w->edges[i], ind2);
    }


    cout << ind << "nb_walls[0]: *\t\t" << (void*)w->nb_walls[0] << " [wall] \t\t/* Array of pointers to walls that share an edge*/ // TODO\n";
    cout << ind << "nb_walls[1]: *\t\t" << (void*)w->nb_walls[1] << " [wall] \t\t/* Array of pointers to walls that share an edge*/ // TODO\n";
    cout << ind << "nb_walls[2]: *\t\t" << (void*)w->nb_walls[2] << " [wall] \t\t/* Array of pointers to walls that share an edge*/ // TODO\n";

    cout << ind << "area: \t\t" << w->area << " [double] \t\t/* Area of this element */\n";

    cout << ind << "normal: \t\t" << w->normal << " [vector3] \t\t/* Normal vector for this wall */\n";
    cout << ind << "unit_u: \t\t" << w->unit_u << " [vector3] \t\t/* U basis vector for this wall */\n";
    cout << ind << "unit_v: \t\t" << w->unit_v << " [vector3] \t\t/* V basis vector for this wall */\n";
    cout << ind << "d: \t\t" << w->d << " [double] \t\t/* Distance to origin (point normal form) */\n";

    cout << ind << "grid: *\t\t" << w->grid << " [surface_grid] \t\t/* Grid of effectors for this wall */\n";

    cout << ind << "flags: \t\t" << w->flags << " [u_short] \t\t/* Count Flags: flags for whether and what we need to count */\n";

    cout << ind << "parent_object: *\t\t" << w->parent_object << " [object] \t\t/* The object we are a part of */\n";
    //dump_object(w->parent_object, ind2);

    cout << ind << "birthplace: *\t\t" << w->birthplace << " [storage] \t\t/* Where we live in memory */\n";

    cout << ind << "counting_regions: *\t\t" << w->counting_regions << " [region_list] \t\t/* Counted-on regions containing this wall */\n";
  }
}


void dump_wall_array(int num, wall** wall_array, const char* ind) {
  if (wall_array == NULL) {
    cout << "\n";
    return;
  }
  for (int i = 0; i < num && i < MAX_ARRAY_ITEMS; i++) {
    if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
      cout << "\n" << ind << "  ";
    }
    cout << ind << i << ":\n";
    dump_wall(wall_array[i], IND_ADD2(ind), false);
  }
}


void dump_wall_list(wall_list* list, const char* ind) {
  wall_list* curr = list;
  int i = 0;
  while (curr != NULL) {
    cout << ind << i << ":\n";
    dump_wall(curr->this_wall, IND_ADD2(ind), false);
    i++;
    curr = curr->next;
  }
}


void dump_wall_list_array(int num, const char* num_name, wall_list** values, const char* values_name, const char* comment, const char* ind) {
  cout << ind << values_name << "[" << num_name << "]: \t\t" << values << "[" << num << "]" << " [wall_list*] \t\t" << comment;
  if (values == NULL) {
    cout << "\n";
    return;
  }
  for (int i = 0; i < num && i < MAX_ARRAY_ITEMS; i++) {
    if (i % DUMP_ARRAY_NEWLINE_COUNT == 0) {
      cout << "\n" << ind << "  ";
    }
    cout << ind << i << ":\n";
    dump_wall_list(values[i], IND_ADD2(ind));
  }
  if (num >= MAX_ARRAY_ITEMS) {
    cout << "...";
  }
  cout << "\n";
}



void dump_molecule_flags(short flags, const char* ind) {
  cout << ind << "flags: \t\t" << flags << " [short] \t\t /* Abstract Molecule Flags: Who am I, what am I doing, etc. */\n";

  cout << ind << "  " << get_molecule_flags_string(flags) << "\n";
}

void dump_abstract_molecule(abstract_molecule* amp, const char* ind) {
  cout << ind << "this: *\t\t" << (void*)amp << " [abstract_molecule] \t\t /* Current molecule */ //???\n";
  cout << ind << "next: *\t\t" << (void*)amp->next << " [abstract_molecule] \t\t /* Next molecule in scheduling queue */ //???\n";
  cout << ind << "t: \t\t" << amp->t << " [double] \t\t                      /* Scheduling time. */\n";
  cout << ind << "t2: \t\t" << amp->t2 << " [double] \t\t                     /* Time of next unimolecular reaction */\n";
  dump_molecule_flags(amp->flags, "  ");

  cout << ind << "properties: *\t\t" << (void*)amp->properties << " [species] \t\t    /* What type of molecule are we? */\n";
  cout << ind << "birthplace: *\t\t" << (void*)amp->birthplace << " [mem_helper] \t\t /* What was I allocated from? */\n";
  cout << ind << "birthday: \t\t" << amp->birthday << " [double] \t\t               /* Real time at which this particle was born */\n";
  cout << ind << "id: \t\t" << amp->id << " [u_long] \t\t                     /* unique identifier of this molecule */\n";
  cout << ind << "periodic_box: *\t\t" << (void*)amp->periodic_box << " [periodic_image] \t\t  /* track the periodic box a molecule is in */\n";

  cout << ind << "graph_data: *\t\t" << (void*)amp->graph_data << " [graph_data] \t\t /* nfsim graph structure data */\n";
#if 0 // nfsim function pointers, not important for now
  u_int (*get_flags)(void *); /* (nfsim) returns the reactivity flags associated with this particle */
  double (*get_diffusion)(void *);        /* (nfsim) returns the diffusion value */
  double (*get_time_step)(void *);        /* (nfsim) function pointer to a method that returns the time step */
  double (*get_space_step)(void *);       /* (nfsim) function pointer to a method that returns the space step */
  /* end structs used by the nfsim integration */
#endif
  cout << ind << "mesh_name: *\t\t" << amp->mesh_name << " [char] \t\t                /* Name of mesh that molecule is either in (volume molecule) or on (surface molecule) */\n";
}


void dump_surface_molecule(surface_molecule* amp, const char* ind) {
  cout << "***TODO!\n";
}

void dump_molecules(int num_all_molecules, molecule_info **all_molecules) {
  cout << "all_molecules: **\t\t" << (void**)all_molecules << " [molecule_info] \n";
  cout << "num_all_molecules: \t\t" << num_all_molecules << " [int] \n";
  for (int i = 0; i < num_all_molecules; i++) {
    molecule_info *curr = all_molecules[i];
    assert(curr != NULL);

    abstract_molecule *amp = curr->molecule;
    assert(amp != NULL);

    if ((amp->properties->flags & NOT_FREE) == 0) {
      volume_molecule *vmp = (volume_molecule *)amp;
      dump_volume_molecule(vmp, "  ", true, "", 0, 0.0, true);
    } else if ((amp->properties->flags & ON_GRID) != 0) {
      surface_molecule *smp = (surface_molecule *)amp;
      dump_surface_molecule(smp, "  ");
    } else {
      dump_abstract_molecule(amp, "  ");
    }

    cout << "reg_names: *\t\t" << curr->reg_names << " [string_buffer] \t\t   /* Region names */\n";
    cout << "mesh_names: *\t\t" << curr->mesh_names << " [string_buffer] \t\t  /* Mesh names that molec is nested in */\n";
    cout << "pos: \t\t" << curr->pos << " [vector3] \t\t /* Position in space */\n";
    cout << "orient: \t\t" << curr->orient << " [short] \t\t /* Which way do we point? */\n";
  }
}


void dump_sm_dat(sm_dat* smd, const char* name, const char* comment, const char* ind) {
  DECL_IND2(ind)
  cout << ind << name << ": *\t\t" << (void*)smd << " [sm_dat] \t\t " << comment << "\n";

  sm_dat* curr = smd;
  int i = 0;
  while (curr != nullptr) {
    cout << ind << i << ":\n";

    cout << ind2 << "sm: *\t\t" << get_sym_name(curr->sm->sym) << " [species name] \t\t /* Species to place on surface */\n";
    cout << ind2 << "quantity_type: *\t\t" << ((curr->quantity_type == SURFMOLDENS) ? "SURFMOLDENS" : "SURFMOLNUM") << " [byte] \t\t  // Placement Type Flags: either SURFMOLDENS or SURFMOLNUM\n";
    cout << ind2 << "quantity: \t\t" << curr->quantity << " [double] \t\t // Amount of surface molecules to place by density or number\n";
    cout << ind2 << "orientation: \t\t" << curr->orientation << " [short] \t\t /* Orientation of molecules to place */\n";

    i++;
    curr = curr->next;
  }
}


void dump_region(region* reg, const char* ind) {
  DECL_IND2(ind)

  // TODO XXXXXXXX
  cout << ind << "reg: *\t\t" << (void*)reg << " [region] \t\t\n";
  if (reg == nullptr) {
    return;
  }
  cout << ind << "  " << "sym: *\t\t" << reg->sym << " [sym_entry] \t\t  /* Symbol hash table entry for this region */\n";
  cout << ind << "  " << "hashval: \t\t" << reg->hashval << " [u_int] \t\t          /* Hash value for counter hash table */\n";
  cout << ind << "  " << "region_last_name: *\t\t" << reg->region_last_name << " [char] \t\t /* Name of region without prepended object name */\n";
  cout << ind << "  " << "parent: *\t\t" << (void*)reg->parent << " [object] \t\t  /* Parent of this region */\n";
  cout << ind << "  " << "element_list_head: *\t\t" << (void*)reg->element_list_head << " [element_list] \t\t /* List of element ranges comprising this region (used at parse time) */\n";
  cout << ind << "  " << "membership: *\t\t" << (void*)reg->membership << " [bit_array] \t\t /* Each bit indicates whether the corresponding wall is in the region */\n";

  //cout << ind << "  " << "sm_dat_head: *\t\t" << (void*)reg->sm_dat_head << " [sm_dat] \t\t /* List of surface molecules to add to region */\n";
  dump_sm_dat(reg->sm_dat_head, "sm_dat_head", "/* List of surface molecules to add to region */", ind2);

  dump_species(reg->surf_class, "surf_class", "/* Surface class of this region */", ind2);
  cout << ind << "  " << "bbox: *\t\t" << (void*)reg->bbox << " [vector3] \t\t /* Array of length 2 to hold corners of region bounding box (used for release in region) */\n";
  cout << ind << "  " << "area: \t\t" << reg->area << " [double] \t\t          /* Area of region */\n";
  cout << ind << "  " << "flags: \t\t" << reg->flags << " [u_short] \t\t        /* Counting subset of Species Flags */\n";
  cout << ind << "  " << "manifold_flag: \t\t" << (int)reg->manifold_flag << " [byte] \t\t   /* Manifold Flags: If IS_MANIFOLD, region is a closed manifold and thus defines a volume */\n";
  cout << ind << "  " << "volume: \t\t" << reg->volume << " [double] \t\t                   /* volume of region for closed manifolds */\n";
  cout << ind << "  " << "boundaries: *\t\t" << (void*)reg->boundaries << " [pointer_hash] \t\t /* hash table of edges that constitute external boundary of the region */\n";
  cout << ind << "  " << "region_has_all_elements: \t\t" << reg->region_has_all_elements << " [int] \t\t /* flag that tells whether the region contains ALL_ELEMENTS (effectively comprises the whole object) */\n";
}

void dump_region_list(region_list *rl, const char* regions_name, const char* comment, const char* ind) {
  cout << ind << regions_name << ": *\t\t" << rl << " [region_list] \t\t " << comment << "\n";
  for (region_list *r = rl; r != NULL; r = r->next) {
    dump_region(r->reg, IND_ADD2(ind) );
  }
}


void dump_one_waypoint(waypoint* wp, const char* ind) {
  assert(wp != NULL);
  cout << ind << "waypoint: *\t\t" << (void*)wp << " [waypoint] \t\t\n";
  cout << ind << "  " << "loc: \t\t" << wp->loc << " [vector3] \t\t          /* This is where the waypoint is */\n";
  dump_region_list(wp->regions, "regions", "/* We are inside these regions */", IND_ADD2(ind));
  dump_region_list(wp->antiregions, "antiregions", "/* We are outside of (but hit) these regions */", IND_ADD2(ind));
}


void dump_waypoints(int n_waypoints, waypoint* waypoints) {
  cout << "n_waypoints: \t\t" << n_waypoints << " [int] \t\t/* How many waypoints (one per subvol) */\n";
  cout << "waypoints: *\t\t" << (void*)waypoints << " [waypoint] \t\t/* Waypoints contain fully-closed region information */\n";

  for (int i = 0; i < n_waypoints; i++) {
    dump_one_waypoint(&waypoints[i], "  ");
  }
}


void dump_one_subvolume(subvolume* sv, const char* ind) {
  assert(sv != NULL);
  DECL_IND2(ind);
  cout << ind << "subvolume: *\t\t" << (void*)sv << " [subvolume] \t\t\n";

#ifdef DUMP_WALLS
  cout << ind << "  " << "wall_head: *\t\t" << (void*)sv->wall_head << " [wall_list] \t\t /* Head of linked list of intersecting walls */\n";
  dump_wall_list(sv->wall_head, ind2);
#endif

  cout << ind << "  " << "mol_by_species: \t\t" << sv->mol_by_species << " [pointer_hash] \t\t /* table of speciesv->molecule list */\n";
  cout << ind << "  " << "species_head: *\t\t" << (void*)sv->species_head << " [per_species_list] \t\t\n";
  cout << ind << "  " << "mol_count: \t\t" << sv->mol_count << " [int] \t\t /* How many molecules are here? */\n";
  cout << ind << "  " << "llf: \t\t" << sv->llf << " [int3D] \t\t /* Indices of left lower front corner */\n";
  cout << ind << "  " << "urb: \t\t" << sv->urb << " [int3D] \t\t /* Indices of upper right back corner */\n";
  cout << ind << "  " << "world_edge: \t\t" << sv->world_edge << " [short] \t\t /* Direction Bit Flags that are set for SSVs at edge of world */\n";
  cout << ind << "  " << "local_storage: *\t\t" << (void*)sv->local_storage << " [storage] \t\t /* Local memory and scheduler */\n";
}


void dump_subvolumes(int n_subvols, const char* num_name, const char* num_comment, subvolume* subvols, const char* subvols_name, const char* name_comment, const char* ind) {
  cout << ind << num_name << ": \t\t" << n_subvols << " [int] \t\t" << num_comment << "\n";
  cout << ind << subvols_name << ": *\t\t" << (void*)subvols << " [subvolume] \t\t" << name_comment << "\n";

  for (int i = 0; i < n_subvols && i < MAX_SUBVOLUMES; i++) {
    dump_one_subvolume(&subvols[i], "  ");
  }
  if (n_subvols >= MAX_SUBVOLUMES) {
    cout << ind << "  " << "... (total " << n_subvols << ")\n";
  }
}


void dump_product_list(product* prod, const char* ind) {
  DECL_IND2(ind);
  product* prod_ptr = prod;
  int i = 0;
  while (prod_ptr != nullptr) {
    cout << ind << i << ":\n";
    cout << ind2 << "orientation: \t\t" << prod_ptr->orientation << " [short] \t\t/* Orientation to place molecule */\n";
    dump_species(prod_ptr->prod, "product", "", ind2);
    prod_ptr = prod_ptr->next;
  }


}

void dump_pathway(pathway* pathway_ptr, const char* ind) {
  assert(pathway_ptr != nullptr);
  //DECL_IND2(ind);

  cout << ind << "next: *\t\t" << (void*)pathway_ptr->next << " [pathway] \t\t/* Next pathway for this reaction */\n";
  cout << ind << "pathname: *\t\t" << (void*)pathway_ptr->pathname << " [rxn_pathname] \t\t/* Data for named reaction pathway or NULL */\n";

  //cout << ind << "reactant1: *\t\t" << pathway_ptr->reactant1 << " [species] \t\t/* First reactant in reaction pathway */\n";
  dump_species(pathway_ptr->reactant1, "reactant1", "/* First reactant in reaction pathway */", ind);
  //cout << ind << "reactant2: *\t\t" << pathway_ptr->reactant2 << " [species] \t\t/* Second reactant (NULL if none) */\n";
  dump_species(pathway_ptr->reactant2, "reactant2", "/* Second reactant in reaction pathway */", ind);
  //cout << ind << "reactant3: *\t\t" << pathway_ptr->reactant3 << " [species] \t\t/* Third reactant (NULL if none) */\n";
  dump_species(pathway_ptr->reactant3, "reactant3", "/* Third reactant in reaction pathway */", ind);

  cout << ind << "km: \t\t" << pathway_ptr->km << " [double] \t\t/* Rate constant */\n";
  cout << ind << "km_filename: *\t\t" << ((pathway_ptr->km_filename == nullptr) ? "NULL" : pathway_ptr->km_filename) << " [char*] \t\t/* Filename for time-varying rates */\n";
  cout << ind << "orientation1: \t\t" << pathway_ptr->orientation1 << " [short] \t\t/* Orientation of first reactant */\n";
  cout << ind << "orientation2: \t\t" << pathway_ptr->orientation2 << " [short] \t\t/* Orientation of second reactant */\n";
  cout << ind << "orientation3: \t\t" << pathway_ptr->orientation3 << " [short] \t\t/* Orientation of third reactant */\n";

  cout << ind << "product_head: *\t\t" << (void*)pathway_ptr->product_head << " [product] \t\t/* Linked lists of species created */\n";
  dump_product_list(pathway_ptr->product_head, ind);


  cout << ind << "prod_signature: *\t\t" << ((pathway_ptr->prod_signature == nullptr) ? "NULL" : pathway_ptr->prod_signature) << " [char*] \t\t/* string created from the names of products put in alphabetical order */\n";
  cout << ind << "flags: \t\t" << pathway_ptr->flags << " [short] \t\t/* flags describing special reactions - REFLECTIVE, TRANSPARENT, CLAMP_CONCENTRATION */\n";
}

void dump_pathway_list(pathway* pathway_head, const char* name,  const char* comment, const char* ind) {
  DECL_IND2(ind);
  if (pathway_head == nullptr) {
    cout << ind << name << ": *\t\t" << (void*)pathway_head << " [pathway*] \t\t" << comment << "\n";
  }
  pathway* pathway_ptr = pathway_head;
  int idx = 0;
  while (pathway_ptr != nullptr) {
    cout << ind << name << "[" << idx << "]" << ": *\t\t" << (void*)pathway_ptr << " [pathway*] \t\t" << comment << "\n";
    dump_pathway(pathway_ptr, ind2);
    pathway_ptr = pathway_ptr->next;
  }
}

void dump_rxn_pathname(rxn_pathname* pathname, const char* ind) {
  DECL_IND2(ind);
  cout << ind << "pathname: *\t\t" << (void*)pathname << " [rxn_pathname]\n";
  if (pathname == nullptr) {
    return;
  }

  cout << ind2 << "sym: *\t\t" << pathname->sym << " [sym_entry] \t\t/* Ptr to symbol table entry for this rxn name */\n";
  cout << ind2 << "hashval: \t\t" << pathname->hashval << " [u_int] \t\t/* Hash value for counting named rxns on regions */\n";
  cout << ind2 << "path_num: \t\t" << pathname->path_num << " [u_int] \t\t/* Pathway number in rxn */\n";
  cout << ind2 << "rx: *\t\t" << (void*)pathname->rx << " [rxn] \t\t/* The rxn associated with this name */\n";
  cout << ind2 << "magic: *\t\t" << (void*)pathname->magic << " [magic_list] \t\t/* A list of stuff that magically happens when the reaction happens */\n";
}

void dump_pathway_infos(int n_pathways, const char* count_name, pathway_info* pathway_info_array, const char* name,  const char* comment, const char* ind) {
  DECL_IND2(ind);
  if (pathway_info_array == nullptr) {
    cout << ind << name << ": *\t\t" << (void*)pathway_info_array << " [pathway*] \t\t" << comment << "\n";
  }

  for (int i = 0; i < n_pathways; i++) {
    cout << ind << name << "[" << i << "]" << ":\n";
    cout << ind2 << "count: *\t\t" << pathway_info_array[i].count << " [double] \t\t/* How many times the pathway has been taken */\n";
    //dump_pathway(&pathway_info_array[n_pathways], ind2);
    dump_rxn_pathname(pathway_info_array[i].pathname, ind2);
  }
}

void dump_variable_rxn_rates(t_func* prob_t, const char* name, const char* ind) {
  DECL_IND2(ind);
  cout << ind << "prob_t: *\t\t" << (void*)prob_t << " [t_func] \t\t/* List of probabilities changing over time, by pathway */\n";

  int i = 0;
  for (t_func* curr = prob_t; curr != nullptr; curr = curr->next) {
    cout << ind2 << i << ": time: " << curr->time << ", value: " << curr->value << ", path: " << curr->path << "\n";
    i++;
  }
}


void dump_molecule_species(abstract_molecule *m) {
  if ((m->properties->flags & EXTERNAL_SPECIES) == 0) {
    cout << get_sym_name(m->properties->sym);
  }
  else {
    assert(m->graph_data != nullptr);
    cout << graph_pattern_to_bngl(m->graph_data->graph_pattern);
  }
}


void dump_rxn_substance(species* s, short orient) {
  if ((s->flags & EXTERNAL_SPECIES) == 0) {
    cout << get_sym_name(s->sym);
  }
  else {
    cout << "BNG";
  }
}


void dump_rxn_pathway_for_diff(pathway* pw) {
  dump_rxn_substance(pw->reactant1, pw->orientation1);
  if (pw->reactant2 != nullptr) {
    cout << " + ";
    dump_rxn_substance(pw->reactant2, pw->orientation2);
  }

  cout << " -> ";

  product* prod = pw->product_head;
  while (prod != nullptr) {
    dump_rxn_substance(prod->prod, prod->orientation);
    prod = prod->next;
  }
  cout << "\n";
}


void dump_rxn_for_diff(rxn* rx) {
  cout << "max_fixed_p: \t\t" << rx->max_fixed_p << " [float_t]\n";
  cout << "min_noreaction_p: \t\t" << rx->min_noreaction_p << " [float_t]\n";

  cout << "cum_probs: ";
  for (int i = 0; i < rx->n_pathways; i++) {
    cout << rx->cum_probs[i] << ", ";
  }
  cout << "\n";

  pathway* current_pathway = rx->pathway_head;
  for (int i = 0; i < rx->n_pathways && current_pathway != nullptr; i++) {
    dump_rxn_pathway_for_diff(current_pathway);
    current_pathway = current_pathway->next;
  }
}

void dump_rxn(rxn* rx, const char* ind, bool for_diff) {

  if (for_diff) {
    dump_rxn_for_diff(rx);
    return;
  }
  cout << ind << "sym: *\t\t" << rx->sym << " [sym_entry] \t\t/* Ptr to symbol table entry for this rxn */\n";

  cout << ind << "n_reactants: \t\t" << rx->n_reactants << " [u_int] \t\t/* How many reactants? (At least 1.) */\n";
  cout << ind << "n_pathways: \t\t" << rx->n_pathways << " [int] \t\t/* How many pathways lead away? (Negative = special reaction, i.e. transparent etc...) */\n";

#if 0
  for (int j = 0; j < rx->n_pathways; ++j) {
    struct rxn_pathname *path = rx->info[j].pathname;

    if (path != NULL && strcmp(path->sym->name, rx_name) == 0) {
      *found_rx = rx;
      // XXX below we convert from u_int to int which is bad
      // unfortunately MCell internally mixes u_int and int for
      // the pathway id.
      *path_id = path->path_num;
      cout << "0: \t\t" << s->0 << " [return] \t\t\n";
    }
  }
#endif


  dump_double_array(rx->n_pathways, "n_pathways", rx->cum_probs, "cum_probs", "/* Cumulative probabilities for (entering) all pathways */", ind);

  cout << ind << "max_fixed_p: \t\t" << rx->max_fixed_p << " [double] \t\t/* Maximum 'p' for region of p-space for all  non-cooperative pathways */\n";
  cout << ind << "min_noreaction_p: \t\t" << rx->min_noreaction_p << " [double] \t\t/* Minimum 'p' for region of p-space which is always in the non-reacting pathway. (note that cooperativity may mean that some values of p less than this still do not produce a reaction) */\n";
  cout << ind << "pb_factor: \t\t" << rx->pb_factor << " [double] \t\t/* Conversion factor from rxn rate to rxn probability (used for cooperativity) */\n";

  cout << ind << "product_idx_aux: *\t\t" << (void*)rx->product_idx_aux << " [int] \t\t/* Number of unique players in each pathway. Used for on-the fly calculation of product_idx indexes TODO: nfsim*/\n";
  //dump_int_array(rx->n_pathways, "n_pathways", rx->product_idx_aux, "product_idx_aux", "/* Number of unique players in each pathway. Used for on-the fly calculation of product_idx indexes */", ind);

  cout << ind << "product_idx: *\t\t" << rx->product_idx << " [u_int] \t\t/* Index of 1st player for products of each pathway TODO: nfsim */\n";
  //dump_int_array(rx->n_pathways, "n_pathways", (int*)rx->product_idx, "product_idx", "/* Index of 1st player for products of each pathway */", ind);

  cout << ind << "players: **\t\t" << (void**)rx->players << " [species] \t\t/* Identities of reactants/products */\n";
  dump_species_list(rx->n_reactants/*might be *n_pathways*/, "n_reactants", rx->players, "players", "/* Identities of reactants/products */", "    ");


  /* this information is kept in a separate array because with nfsim we dont know ahead of time how many paths/products per path are there, we only know when it fires*/
  cout << ind << "nfsim_players: ***\t\t" << (void***)rx->nfsim_players << " [species] \t\t/* a matrix of the nfsim elements associated with each path */\n";
  cout << ind << "geometries: *\t\t" << (void*)rx->geometries << " [short] \t\t/* Geometries of reactants/products */\n";
  cout << ind << "nfsim_geometries: **\t\t" << (void**)rx->nfsim_geometries << " [short] \t\t/* geometries of the nfsim geometries associated with each path */\n";

  cout << ind << "n_occurred: \t\t" << rx->n_occurred << " [long long] \t\t/* How many times has this reaction occurred? */\n";
  cout << ind << "n_skipped: \t\t" << rx->n_skipped << " [double] \t\t/* How many reactions were skipped due to probability overflow? */\n";

  dump_variable_rxn_rates(rx->prob_t, "prob_t", "    ");


  cout << ind << "pathway_head: *\t\t" << rx->pathway_head << " [pathway] \t\t/* List of pathways built at parse-time */\n";

  dump_pathway_list(rx->pathway_head, "pathway_head", "/* List of pathways built at parse-time */", "    ");

  cout << ind << "info: *\t\t" << (void*)rx->info << " [pathway_info] \t\t/* Counts and names for each pathway */\n";

  dump_pathway_infos(rx->n_pathways, "n_pathways", rx->info, "info", "/* Counts and names for each pathway */", "    ");

  cout << ind << "TODO: NFSIM callbacks\n";
#if 0
  //char** external_reaction_names; /* Stores reaction results stored from an external program (like nfsim)*/
  external_reaction_datastruct* external_reaction_data; /* Stores reaction results stored from an external program (like nfsim)*/
  graph_data** reactant_graph_data; /* stores the graph patterns associated with the reactants for every path */
  graph_data*** product_graph_data; /* Stores the graph patterns associated with our products for each path*/
  double (*get_reactant_diffusion)(rxn*, int);  /* returns the diffusion value associated with its reactants*/
  double (*get_reactant_time_step)(rxn*, int);  /* returns the diffusion value associated with its reactants*/
  double (*get_reactant_space_step)(rxn*, int);  /* returns the diffusion value associated with its reactants*/
#endif


}

void dump_reaction_hash_table(int rx_hashsize, const char* num_name, rxn **reaction_hash, const char* reaction_hash_name, const char* comment, const char* ind) {
  cout << ind << reaction_hash_name << "[" << num_name << "]: \t\t" << (void*)reaction_hash << "[" << rx_hashsize << "]" << " [reaction_hash*] \t\t" << comment << "\n";
  // used code from get_rxn_by_name
  for (int i = 0; i < rx_hashsize; ++i) {
    struct rxn *rx = NULL;
    struct rxn *rx_array = reaction_hash[i];
    cout << IND_ADD2(ind) << i << ": " << (void*)rx_array << "\n";
    int k = 0;
    for (rx = rx_array; rx != NULL; rx = rx->next) { // go over linked list
      cout << IND_ADD2(IND_ADD2(ind)) << k << ":\n";
      dump_rxn(rx, IND_ADD2(IND_ADD2(ind)));
    }
  }
}


void dump_region_data_owners(
    int n_objects, int* obj_index, geom_object** owners, const char* name, const char* comments, const char* ind
) {
  cout << ind << "owners: **\t\t" << owners << " [object] \t\t/* Array of pointers to each object */\n";

  if (owners == nullptr) {
    return;
  }
  DECL_IND2(ind);

  // get max index for owners array
  int max_idx = -1;
  for (int i = 0; i < n_objects; i++) {
    if (obj_index[i] > max_idx) {
      max_idx = obj_index[i];
    }
  }


  for (int i = 0; i <= max_idx; i++) {
    cout << ind2 << i << ":";
    dump_object(owners[i], ind2);
  }
}


// recursive
void dump_release_evaluator(release_evaluator* expr, const char* ind) {
  DECL_IND2(ind);
  if (expr == nullptr) {
    return;
  }

  cout << ind2 << "op: " << get_region_expression_flags_string(expr->op) << "\n";
  cout << ind2 << "left: ";
  if (expr->op & REXP_LEFT_REGION) {
    region* r = (region *)expr->left;
    dump_region(r, ind2);
  }
  else {
    dump_release_evaluator((release_evaluator*)expr->left, ind2);
  }
  cout << ind2 << "right: ";
  if (expr->op & REXP_RIGHT_REGION) {
    region* r = (region *)expr->right;
    dump_region(r, ind2);
  }
  else {
    dump_release_evaluator((release_evaluator*)expr->right, ind2);
  }
}


void dump_release_region_data(release_region_data* region_data, const char* name, const char* comment, const char* ind) {
  DECL_IND2(ind);

  cout << ind << name << ": *\t\t" << (void*)region_data << " [release_region_data] \t\t" << comment << "\n";
  if (region_data != nullptr) {
    cout << ind2 << "llf: \t\t" << region_data->llf << " [vector3] \t\t/* One corner of bounding box for release volume */\n";
    cout << ind2 << "urb: \t\t" << region_data->urb << " [vector3] \t\t/* Opposite corner */\n";

    //cout << ind2 << "n_walls_included: \t\t" << region_data->n_walls_included << " [int] \t\t/* How many walls total */\n";
    //cout << ind2 << "cum_area_list: *\t\t" << region_data->cum_area_list << " [double] \t\t/* Cumulative area of all walls */\n";
    dump_double_array(region_data->n_walls_included, "n_walls_included", region_data->cum_area_list, "cum_area_list", "/* Cumulative area of all walls */", ind2);

    //cout << ind2 << "wall_index: *\t\t" << region_data->wall_index << " [int] \t\t/* Indices of each wall (by object) */\n";
    dump_int_array(region_data->n_walls_included, "n_walls_included", region_data->wall_index, "wall_index", "/* Indices of each wall (by object) */", ind2);

    cout << ind2 << "n_objects: \t\t" << region_data->n_objects << " [int] \t\t/* How many objects are there total */\n";
    //cout << ind2 << "obj_index: *\t\t" << region_data->obj_index << " [int] \t\t/* Indices for objects (in owners array) */\n";
    dump_int_array(region_data->n_walls_included, "n_walls_included", region_data->obj_index, "obj_index", "/* count: n_objects Indices for objects (in owners array) */", ind2);

    // TODO - dump
    dump_region_data_owners(region_data->n_objects, region_data->obj_index, region_data->owners, "owners", "Array of pointers to each object", ind2);

    // TODO - dump
    cout << ind2 << "in_release: **\t\t" << region_data->in_release << " [bit_array] \t\t/* Array of bit arrays; each bit array says which walls are in release for an object */\n";
    // TODO - dump
    cout << ind2 << "walls_per_obj: *\t\t" << region_data->walls_per_obj << " [int] \t\t/* Number of walls in release for each object */\n";

    cout << ind2 << "self: *\t\t" << region_data->self << " [object] \t\t/* A pointer to our own release site object */\n";
    dump_object(region_data->self,  IND_ADD2(ind2));

    cout << ind2 << "expression: *\t\t" << region_data->expression << " [release_evaluator] \t\t/* A set-construction expression combining regions to form this release site */\n";
    dump_release_evaluator(region_data->expression, ind2);
  }

}


void dump_release_pattern(release_pattern* pattern, const char* name, const char* comment, const char* ind) {
  DECL_IND2(ind);

  cout << ind << name << ": *\t\t" << (void*)pattern << " [release_pattern] \t\t" << comment << "\n";
  if (pattern != nullptr) {
    //cout << ind2 << "sym: *\t\t" << pattern->sym << " [sym_entry] \t\t/* Symbol hash table entry for the pattern */\n";
    cout << ind2 << "sym.name: \t\t" << get_sym_name(pattern->sym) << " [sym_entry.name] \t\t/* Symbol hash table entry for the pattern */\n";
    cout << ind2 << "delay: \t\t" << pattern->delay << " [double] \t\t/* Delay between time 0 and first release event. */\n";
    cout << ind2 << "release_interval: \t\t" << pattern->release_interval << " [double] \t\t/* Time between release events within a train. */\n";
    cout << ind2 << "train_interval: \t\t" << pattern->train_interval << " [double] \t\t/* Time from the start of one train to the start of the next one. */\n";
    cout << ind2 << "train_duration: \t\t" << pattern->train_duration << " [double] \t\t/* Length of the train. */\n";
    cout << ind2 << "number_of_trains: \t\t" << pattern->number_of_trains << " [int] \t\t/* How many trains are produced. */\n";
  }
}


void dump_release_site_obj(release_site_obj* rel_site, const char* ind) {
  DECL_IND2(ind);

  if (rel_site->location != nullptr) {
    cout << ind2 << "location: *\t\t" << *rel_site->location << " [vector3] \t\t/* location of release site */\n";
  }
  else {
    cout << ind2 << "location: *\t\t" << (void*)rel_site->location << " [vector3*] \t\t/* location of release site */\n";
  }
  cout << ind2 << "mol_type: *\t\t" << (void*)rel_site->mol_type << " [species] \t\t/* species to be released */\n";
  //dump_species(rel_site->mol_type, "mol_type", "/* species to be released */", ind2);
  if (rel_site->mol_type != nullptr) {
    cout << ind2 << "  mol_type (name): *\t\t" << rel_site->mol_type->sym->name << "\n";
  }
  else {
    cout << ind2 << "  mol_type: *\t\t" << "NULL " << "\n";
  }

  cout << ind2 << "release_number_method: \t\t" << get_release_number_method_string((release_number_type_t)rel_site->release_number_method) << " [byte] \t\t/* Release Number Flags: controls how release_number is used (enum release_number_type_t) */\n";
  cout << ind2 << "release_shape: \t\t" << get_release_shape_name(rel_site->release_shape) << " [int8_t] \t\t/* Release Shape Flags: controls shape over which to release (enum release_shape_t) */\n";
  cout << ind2 << "orientation: \t\t" << rel_site->orientation << " [short] \t\t/* Orientation of released surface molecules */\n";
  cout << ind2 << "release_number: \t\t" << rel_site->release_number << " [double] \t\t/* Number to release */\n";
  cout << ind2 << "mean_diameter: \t\t" << rel_site->mean_diameter << " [double] \t\t/* Diameter for symmetric releases */\n";
  cout << ind2 << "concentration: \t\t" << rel_site->concentration << " [double] \t\t/* Concentration of molecules to release. Units are Molar for volume molecules, and number per um^2 for surface molecules. */\n";
  cout << ind2 << "standard_deviation: \t\t" << rel_site->standard_deviation << " [double] \t\t/* Standard deviation of release_number for GAUSSNUM, or of mean_diameter for VOLNUM */\n";
  if (rel_site->diameter != nullptr) {
    cout << ind2 << "diameter: *\t\t" << *rel_site->diameter << " [vector3] \t\t/* x,y,z diameter for geometrical release shapes */\n";
  }
  else {
    cout << ind2 << "diameter: *\t\t" << (void*)rel_site->diameter << " [vector3*] \t\t/* x,y,z diameter for geometrical release shapes */\n";
  }

#ifdef DUMP_RELEASE_REGION_DATA
  dump_release_region_data(rel_site->region_data, "region_data", "/* Information related to release on regions */", ind2);
#else
  cout << ind2 << "region_data: *\t\t" << (void*)rel_site->region_data << " [release_region_data] \t\t/* Information related to release on regions */\n";
#endif

  cout << ind2 << "mol_list: *\t\t" << (void*)rel_site->mol_list << " [release_single_molecule] \t\t/* Information related to release by list */\n";
  cout << ind2 << "release_prob: \t\t" << rel_site->release_prob << " [double] \t\t/* Probability of releasing at scheduled time */\n";

  if (rel_site->periodic_box != nullptr) {
    cout << ind2 << "periodic_box: *\t\t" << *rel_site->periodic_box << " [periodic_image] \t\t\n";
  }
  else {
    cout << ind2 << "periodic_box: *\t\t" << (void*)rel_site->periodic_box << " [periodic_image] \t\t\n";
  }

  dump_release_pattern(rel_site->pattern, "pattern", "/* Timing of releases by virtual function generator */", ind2);

  cout << ind2 << "name: *\t\t" << rel_site->name << " [char] \t\t/* Fully referenced name of the instantiated release_site */\n";
  cout << ind2 << "graph_pattern: *\t\t" << (void*)rel_site->graph_pattern << " [char] \t\t/* JJT: Graph definition of the structured molecules being released in this site*/\n";
}

void dump_release_event_queue(release_event_queue* req, const char* ind) {

  DECL_IND2(ind);
  cout << ind << "release_event_queue :" << (void*)req << "\n";
  cout << ind2 << "event_time: \t\t" << req->event_time << " [double] \t\t/* Time of the release */\n";
  cout << ind2 << "release_site: *\t\t" << (void*)req->release_site << " [release_site_obj] \t\t/* What to release, where to release it, etc */\n";
  dump_release_site_obj(req->release_site, ind2);
  cout << ind2 << "t_matrix: *\t\t" << req->t_matrix << " [double[4][4]] \t\t // transformation matrix for location of release site\n";
  cout << ind2 << "train_counter: \t\t" << req->train_counter << " [int] \t\t/* counts executed trains */\n";
  cout << ind2 << "train_high_time: \t\t" << req->train_high_time << " [double] \t\t/* time of the train's start */\n";

  cout << ind2 << "next: \t\t" << (void*)req->next << "\n";

  //if (req->next != NULL)
  //  dump_release_event_queue(req->next, ind);
}

void dump_schedule_helper(schedule_helper* shp, const char* name, const char* comment, const char* ind, bool simplified_for_vm) {
  if (simplified_for_vm) {
    cout << "Scheduler '" << name << "', now: " << shp->now << "\n    ";

    for (schedule_helper* shp_curr = shp; shp_curr != NULL; shp_curr = shp_curr->next_scale) {

      if (shp_curr->next_scale != NULL) {
        cout << "\nnext_scale:\n";
      }

      abstract_molecule* am = (abstract_molecule*)shp_curr->current;

      while (am != NULL) {
        if (am->properties == NULL) {
          cout << "!NULL properties!,";
        }
        else {
          volume_molecule* vm = (volume_molecule*)am;
          cout << "(t: " << vm->t << ", t2: " << vm->t2 << ", id: " << vm->id << "), ";
        }
        am = am->next;
      }
      cout << "\n";
    }
  }
  else {
    std::string inds = ind;
    inds += "  ";
    const char* ind2 = inds.c_str();
    cout << ind << name << ": *\t\t" << (void*)shp << " [schedule_helper] \t\t" << comment << "\n";
    if (shp == nullptr) {
      return;
    }
#ifdef DUMP_SCHEDULERS
    cout << ind2 <<"next_scale: *\t\t" << (void*)shp->next_scale << " [schedule_helper] \t\t/* Next coarser time scale */\n";

    cout << ind2 <<"dt: \t\t" << shp->dt << " [double] \t\t/* Timestep per slot */\n";
    cout << ind2 <<"dt_1: \t\t" << shp->dt_1 << " [double] \t\t/* dt_1 = 1/dt */\n";
    cout << ind2 <<"now: \t\t" << shp->now << " [double] \t\t/* Start time of the scheduler */\n";

    /* Items scheduled now or after now */
    cout << ind2 <<"count: \t\t" << shp->count << " [int] \t\t/* Total number of items scheduled now or after */\n";
    cout << ind2 <<"buf_len: \t\t" << shp->buf_len << " [int] \t\t/* Number of slots in the scheduler */\n";
    cout << ind2 <<"index: \t\t" << shp->index << " [int] \t\t/* index of the next time block */\n";

    if (shp->circ_buf_count != NULL) {
      cout << ind2 <<"circ_buf_count: \t\t" << *shp->circ_buf_count << " [int] \t\t/* How many items are scheduled in each slot */\n";
    }
    else {
      cout << ind2 <<"circ_buf_count: \t\t" << "NULL" << " [int*] \t\t/* How many items are scheduled in each slot */\n";
    }

    cout << ind2 <<"circ_buf_head: **\t\t" << (void**)shp->circ_buf_head << " [abstract_element] \t\t// Array of linked lists of scheduled items for each slot\n";
    cout << ind2 <<"circ_buf_tail: **\t\t" << (void**)shp->circ_buf_tail << " [abstract_element] \t\t// Array of tails of the linked lists\n";

    cout << ind2 <<"contents (current, circ_buf_head):\n";
    for (int i = -1; i < shp->buf_len; i++) {
      int k = 0;
      for (struct abstract_element *aep = (i < 0) ? shp->current
                                                  : shp->circ_buf_head[i];
           aep != NULL; aep = aep->next) {

        cout << ind2 << "  " << i << ":\n";
        if (strcmp(name, "releaser") == 0) {
          struct release_event_queue *req = (struct release_event_queue *)aep;
          dump_release_event_queue(req, ind2);

        }
        else if (strcmp(name, "dynamic_geometry_scheduler") == 0) {
          // TODO
        }
        else {
          struct abstract_molecule *amp = (struct abstract_molecule *)aep;
          if (amp->properties == NULL) {
            cout << ind2 << "  " << i << "." << k << ": " << (void*)amp << ", properties: " << (void*)amp->properties << "\n";
            k++;
            continue;
          }
          else {
            k++;
          }

          dump_abstract_molecule(amp, IND_ADD2(ind2));
        }
      }
    }

    /* Items scheduled before now */
    /* These events must be serviced before simulation can advance to now */
    cout << ind2 <<"current_count: \t\t" << shp->current_count << " [int] \t\t/* Number of current items */\n";
    cout << ind2 <<"current: *\t\t" << (void*)shp->current << " [abstract_element] \t\t/* List of items scheduled now */\n";
    cout << ind2 <<"current_tail: *\t\t" << shp->current_tail << " [abstract_element] \t\t/* Tail of list of items */\n";

    cout << ind2 <<"defunct_count: \t\t" << shp->defunct_count << " [int] \t\t/* Number of defunct items (set by user)*/\n";
    cout << ind2 <<"error: \t\t" << shp->error << " [int] \t\t/* Error code (1 - on error, 0 - no errors) */\n";
    cout << ind2 <<"depth: \t\t" << shp->depth << " [int] \t\t/* Tier of scheduler in timescale hierarchy, 0-based */\n";
#endif
  }
}



void dump_species_item(species* spec, const char* ind) {
  DECL_IND2(ind);

  cout << ind2 <<"species_id: \t\t" << spec->species_id << " [u_int] \t\t/* Unique ID for this species */\n";
  cout << ind2 <<"chkpt_species_id: \t\t" << spec->chkpt_species_id << " [u_int] \t\t/* Unique ID for this species from the checkpoint file */\n";
  cout << ind2 <<"hashval: \t\t" << spec->hashval << " [u_int] \t\t/* Hash value (may be nonunique) */\n";
  cout << ind2 <<"sym: *\t\t" << spec->sym << " [sym_entry] \t\t/* Symbol table entry (name) */\n";
  cout << ind2 <<"sm_dat_head: *\t\t" << spec->sm_dat_head << " [sm_dat] \t\t/* If IS_SURFACE this points to head of effector data list associated with surface class */\n";

  cout << ind2 <<"population: \t\t" << spec->population << " [u_int] \t\t/* How many of this species exist? */\n";

  cout << ind2 <<"D: \t\t" << spec->D << " [double] \t\t/* Diffusion constant */\n";
  cout << ind2 <<"space_step: \t\t" << spec->space_step << " [double] \t\t/* Characteristic step length */\n";
  cout << ind2 <<"time_step: \t\t" << spec->time_step << " [double] \t\t/* Minimum (maximum?) sensible timestep */\n";
  cout << ind2 <<"max_step_length: \t\t" << spec->max_step_length << " [double] \t\t/* maximum allowed random walk step */\n";

  cout << ind2 <<"flags: \t\t" << spec->flags << " [u_int] " << get_species_flags_string(spec->flags) << " \t\t/* Species Flags:  Vol Molecule? Surface Molecule? Surface Class? Counting stuff, etc... */\n";

  cout << ind2 <<"n_deceased: \t\t" << spec->n_deceased << " [long long] \t\t/* Total number that have been destroyed. */\n";
  cout << ind2 <<"cum_lifetime_seconds: \t\t" << spec->cum_lifetime_seconds << " [double] \t\t/* Seconds lived by now-destroyed molecules */\n";

  /* if species s a surface_class (IS_SURFACE) below there are linked lists of
   * molecule names/orientations that may be present in special reactions for
   * this surface class */
  cout << ind2 <<"refl_mols: *\t\t" << spec->refl_mols << " [name_orient] \t\t// names of the mols that REFLECT from surface\n";
  cout << ind2 <<"transp_mols: *\t\t" << spec->transp_mols << " [name_orient] \t\t/* names of the mols that are TRANSPARENT for surface */\n";
  cout << ind2 <<"absorb_mols: *\t\t" << spec->absorb_mols << " [name_orient] \t\t// names of the mols that ABSORB at surface\n";
  cout << ind2 <<"clamp_conc_mols: *\t\t" << spec->clamp_conc_mols << " [name_orient] \t\t/* names of mols that CLAMP_CONC at surface */\n";
}


void dump_species(species* spec, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << (void*)spec << " [species] \t\t" << comment << "\n";
  if (spec != nullptr) {
    dump_species_item(spec, ind);
  }
}



void dump_species_list(int n_species, const char* num_name, species** species_list, const char* name, const char* comment, const char* ind) {
  DECL_IND2(ind);

  cout << ind2 << name << "[" << num_name << "]: \t\t" << (void**)species_list << "[" << n_species << "]" << " [species*[]] \t\t" << comment << "\n";

  for (int i = 0; i < n_species; i++) {
    species* spec = species_list[i];
    cout << ind2 << "  " << i << ": " << (void*)spec << "\n";
    dump_species_item(spec, IND_ADD2(ind2));
  }
}



string viz_output_flag_to_str(u_short flag) {
  string res = "";

  if (flag & VIZ_ALL_MOLECULES) {
    res += "VIZ_ALL_MOLECULES";
  }
  if (flag & VIZ_MOLECULES_STATES) {
    res += "VIZ_MOLECULES_STATES";
  }
  if (flag & VIZ_SURFACE_STATES) {
    res += "VIZ_SURFACE_STATES";
  }

  if (res == "") {
    res = "NONE";
  }
  return res;
}



void dump_frame_data_list(frame_data_list* frame_data_head, const char* name, const char* comment, const char* ind) {
  DECL_IND2(ind);
  cout << ind << name << ": *\t\t" << frame_data_head << " [frame_data_list] \t\t" << comment << "\n";

  if (frame_data_head == nullptr) {
    cout << ind2 << "NULL\n";
    return;
  }

  cout << ind2 << "next: *\t\t" << frame_data_head->next << " [frame_data_list] \t\t\n";

  cout << ind2 << "list_type: \t\t" << frame_data_head->list_type << " [output_timer_type_t] \t\t/* Data Output Timing Type (OUTPUT_BY_TIME_LIST, etc) */\n";
  cout << ind2 << "type: \t\t" << frame_data_head->type << " [viz_frame_type_t] \t\t/* Visualization Frame Data Type (ALL_FRAME_DATA, etc) */\n";
  cout << ind2 << "viz_iteration: \t\t" << frame_data_head->viz_iteration << " [long long] \t\t/* Value of the current iteration step. */\n";
  cout << ind2 << "n_viz_iterations: \t\t" << frame_data_head->n_viz_iterations << " [long long] \t\t/* Number of iterations in the iteration_list. */\n";

  cout << ind2 << "iteration_list: *\t\t" << frame_data_head->iteration_list << " [num_expr_list] \t\t/* Linked list of iteration steps values */\n";
  cout << ind2 << "curr_viz_iteration: *\t\t" << frame_data_head->curr_viz_iteration << " [num_expr_list] \t\t/* Points to the current iteration in the linked list */\n";
}

void dump_viz_output_block(viz_output_block* viz_blocks, const char* name, const char* comment, const char* ind) {
  DECL_IND2(ind);
  cout << ind << name << ": *\t\t" << (void*)viz_blocks << " [viz_output_block] \t\t" << comment << "\n";
  if (viz_blocks == nullptr) {
    return;
  }

  cout << ind2 << "next: *\t\t" << viz_blocks->next << " [viz_output_block] \t\t/* Link to next block */\n";

  dump_frame_data_list(viz_blocks->frame_data_head, "frame_data_head", "/* head of the linked list of viz frames to output */", ind2);


  cout << ind2 << "viz_mode: \t\t" << viz_blocks->viz_mode << " [viz_mode_t] \t\t\n";
  cout << ind2 << "file_prefix_name: *\t\t" << viz_blocks->file_prefix_name << " [char] \t\t\n";

  cout << ind2 << "viz_output_flag: \t\t" << viz_output_flag_to_str(viz_blocks->viz_output_flag) << " [u_short] \t\t/* Takes VIZ_ALL_MOLECULES, VIZ_MOLECULES_STATES, etc. */\n";


  if (viz_blocks->species_viz_states == nullptr) {
    cout << ind2 << "species_viz_states: *\t\t" << viz_blocks->species_viz_states << " [int (ptr)] \t\t\n";
  }
  else {
    cout << ind2 << "species_viz_states: \t\t" << *viz_blocks->species_viz_states << " [int] \t\t\n";
  }

  cout << ind2 << "default_mol_state: \t\t" << viz_blocks->default_mol_state << " [int] \t\t// Only set if (viz_output_flag & VIZ_ALL_MOLECULES)\n";

  /* Parse-time only: Tables to hold temporary information. */
  //struct pointer_hash parser_species_viz_states;

}


void dump_output_trigger_data(output_trigger_data* otd, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << otd << " [output_trigger_data] \t\t" << comment << "\n";
  if (otd == nullptr) {
    return;
  }
  DECL_IND2(ind);

  cout << ind2 << "t_iteration: \t\t" << otd->t_iteration << " [double] \t\t /* Iteration time of the triggering event (in sec) */\n";
  cout << ind2 << "event_time: \t\t" << otd->event_time << " [double] \t\t  /* Exact time of the  event */\n";
  cout << ind2 << "loc: \t\t" << otd->loc << " [vector3] \t\t /* Position of event */\n";
  cout << ind2 << "how_many: \t\t" << otd->how_many << " [int] \t\t       /* Number of events */\n";
  cout << ind2 << "orient: \t\t" << otd->orient << " [short] \t\t       /* Orientation information */\n";
  // TODO: dump flags
  cout << ind2 << "flags: \t\t" << otd->flags << " [short] \t\t        /* Output Trigger Flags */\n";
  if (otd->name != nullptr) {
    cout << ind2 << "name: \t\t" << otd->name << " [char*] \t\t         /* Name to give event */\n";
  }
  else {
    cout << ind2 << "name: \t\t" << (void*)otd->name << " [char*] \t\t         /* Name to give event */\n";
  }
  cout << ind2 << "id: \t\t" << otd->id << " [u_long] \t\t /**/\n";

}


void dump_output_buffer(output_buffer* obuf, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << obuf << " [output_buffer] \t\t" << comment << "\n";
  if (obuf == nullptr) {
    return;
  }
  DECL_IND2(ind);

  cout << ind2 << "data_type: \t\t" << obuf->data_type << " [count_type_t] \t\t /**/\n";

  switch(obuf->data_type) {
  case COUNT_UNSET: // probably wrong
    cout << ind2 << "cval: \t\t" << obuf->val.cval << " [char] \t\t  /**/\n";
    break;
  case COUNT_DBL:
    cout << ind2 << "dval: \t\t" << obuf->val.dval << " [double] \t\t  /**/\n";
    break;
  case COUNT_INT:
    cout << ind2 << "ival: \t\t" << obuf->val.ival << " [int] \t\t  /**/\n";
    break;
  case COUNT_TRIG_STRUCT:
    //cout << ind2 << "tval: \t\t" << obuf->val.tval << " [output_trigger_data*] \t\t  /**/\n";
    dump_output_trigger_data(obuf->val.tval, "val.tval", "", ind2);
    break;
  default:
    assert(false);
  }

}

void dump_one_output_column(output_column* column, const char* ind) {
  cout << ind << "next: \t\t" << (void*)column->next << " [output_column*] \t\t  /* Next column in this set */\n";
  cout << ind << "set: \t\t" << (void*)column->set << " [output_set*] \t\t      /* Which set do we belong to? */\n";
  cout << ind << "initial_value: \t\t" << column->initial_value << " [double] \t\t        /* To continue existing cumulative counts--not implemented yet--and keep track of triggered data */\n";

  //cout << ind << "buffer: \t\t" << column->buffer << " [output_buffer*] \t\t /* Output buffer array (cast based on data_type) */\n";
  dump_output_buffer(column->buffer, "buffer", "/* Data for one output column *", ind);

  cout << ind << "expr: \t\t" << column->expr << " (" << (void*)column->expr << ") [output_expression*] \t\t /* Evaluate this to calculate our value (NULL if trigger) */\n";
}

void dump_output_column(output_column* ocol, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << ocol << " [output_column] \t\t" << comment << "\n";

  DECL_IND2(ind);
  output_column* curr = ocol;
  while (curr != nullptr) {
    dump_one_output_column(curr, ind2);
    curr = curr->next;
  }
}


void dump_one_output_set(output_set* block, const char* ind) {
  cout << ind << "next: \t\t" << (void*)block->next << " [output_set*] \t\t            /* Next data set in this block */\n";
  cout << ind << "block: \t\t" << (void*)block->block << " [output_block*] \t\t         /* Which block do we belong to? */\n";
  cout << ind << "outfile_name: \t\t" << block->outfile_name << " [char*] \t\t                 /* Filename */\n";
  cout << ind << "file_flags: \t\t" << block->file_flags << " [overwrite_policy_t] \t\t /* Overwrite Policy Flags: tells us how to handle existing files */\n";
  cout << ind << "chunk_count: \t\t" << block->chunk_count << " [u_int] \t\t    /* Number of buffered output chunks processed */\n";

  if (block->header_comment != nullptr) {
    cout << ind << "header_comment: \t\t" << block->header_comment << " [char*] \t\t /* Comment character(s) for header */\n";
  }
  else {
    cout << ind << "header_comment: \t\t" << (void*)block->header_comment << " [char*] \t\t /* Comment character(s) for header */\n";
  }
  cout << ind << "exact_time_flag: \t\t" << block->exact_time_flag << " [int] \t\t  /* Boolean value; nonzero means print exact time in TRIGGER statements */\n";

  //cout << ind << "column_head: \t\t" << block->column_head << " [output_column*] \t\t /* Data for one output column */\n";
  dump_output_column(block->column_head, "column_head", "/* Data for one output column */", ind);
}


void dump_output_set(output_set* oset, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << oset << " [output_set] \t\t" << comment << "\n";

  DECL_IND2(ind);
  output_set* curr = oset;

  int i = 0;
  while (curr != nullptr) {
    cout << ind << i << ": \n";
    i++;
    dump_one_output_set(curr, ind2);
    curr = curr->next;
  }
}

void dump_one_output_block(output_block* block, const char* ind) {

  cout << ind << "next: \t\t" << (void*)block->next << " [output_block*] \t\t /* Next in world or scheduler */\n";
  cout << ind << "t: \t\t" << block->t << " [double] \t\t                  /* Scheduled time to update counters */\n";

  cout << ind << "timer_type: \t\t" << block->timer_type << " [output_timer_type_t] \t\t /* Data Output Timing Type (OUTPUT_BY_STEP, etc) */\n";

  cout << ind << "step_time: \t\t" << block->step_time << " [double] \t\t     /* Output interval (seconds) */\n";

  cout << ind << "time_list_head: \t\t" << block->time_list_head << " [num_expr_list*] \t\t /* List of output times/iteration numbers */\n";
  cout << ind << "time_now: \t\t" << block->time_now << " [num_expr_list*] \t\t /* Current entry in list */\n";

  cout << ind << "buffersize: \t\t" << block->buffersize << " [u_int] \t\t   /* Size of output buffer */\n";
  cout << ind << "trig_bufsize: \t\t" << block->trig_bufsize << " [u_int] \t\t /* Size of output buffer for triggers */\n";
  cout << ind << "buf_index: \t\t" << block->buf_index << " [u_int] \t\t    /* Index into buffer (for non-triggers) */\n";

  dump_double_array(block->buffersize, "buffersize", block->time_array, "time_array", "/* Array of output times (for non-triggers) */", ind);

  dump_output_set(block->data_set_head, "data_set_head", "/* Linked list of data sets (separate files) */", ind);
}


void dump_output_blocks(output_block* output_block_head, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << output_block_head << " [output_block] \t\t" << comment << "\n";

  DECL_IND2(ind);
  output_block* curr = output_block_head;

  int i = 0;
  while (curr != nullptr) {
    cout << ind << i << ": \n";
    i++;
    dump_one_output_block(curr, ind2);
    curr = curr->next;
  }
}



void dump_one_output_request(output_request* req, const char* ind) {
  cout << ind << "next: \t\t" << (void*)req->next << " [output_request*] \t\t         /* Next request in global list */\n";
  cout << ind << "requester: \t\t" << req->requester << " (" << (void*)req->requester << ") [output_expression*] \t\t /* Expression in which we appear */\n";
  cout << ind << "count_target: \t\t" << req->count_target << " [sym_entry*] \t\t      /* Mol/rxn we're supposed to count */\n";
  cout << ind << "count_orientation: \t\t" << req->count_orientation << " [short] \t\t             /* orientation of the molecule we are supposed to count */\n";
  cout << ind << "count_location: \t\t" << req->count_location << " [sym_entry*] \t\t    /* Object or region on which we're supposed to count it */\n";
  cout << ind << "report_type: \t\t" << get_report_type_flags_string(req->report_type) << " [byte] \t\t                    /* Output Report Flags telling us how to count */\n";
  cout << ind << "periodic_box: \t\t" << (void*)req->periodic_box << " [periodic_image*] \t\t /* periodic box we are counting in; NULL means that we don't care and count everywhere */\n";
}


void dump_output_requests(output_request* output_request_head, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << output_request_head << " [output_request] \t\t" << comment << "\n";

  DECL_IND2(ind);
  output_request* curr = output_request_head;
  int i = 0;
  while (curr != nullptr) {
    cout << ind << i << ": \n";
    i++;
    dump_one_output_request(curr, ind2);
    curr = curr->next;
  }
}


void dump_sym_table(sym_table_head* t, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << t << " [sym_table_head] \t\t" << comment << "\n";
  if (t == nullptr) {
    return;
  }

  DECL_IND2(ind);
  cout << ind2 << "n_bins: \t\t" << t->n_bins << " [int] \t\t /**/\n";
  cout << ind2 << "n_entries: \t\t" << t->n_entries << " [int] \t\t /**/\n";

  int dumped_symbols = 0;
  for (int i = 0; i < t->n_bins; i++) {

    if (t->entries[i] != nullptr) {
      cout << ind2 << "[" << i << "]:" <<  t->entries[i] << "\n";
      dumped_symbols++;
    }
  }

  if (dumped_symbols != t->n_entries) {
    cout << "internal warning: when dumping symbol table: dumped_symbols != t->n_entries\n";
  }
}


void dump_name_list(name_list* obj, const char* name, const char* comment, const char* ind) {
  name_list* curr = obj;
  int i = 0;
  while (curr != NULL) {
    cout << ind << i << ": " << str(curr->name) << "\n";
    i++;
    curr = curr->next;
  }
}


void dump_dyngeom_parse_vars(dyngeom_parse_vars* dg_parse, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << (void*)dg_parse << " [dyngeom_parse_vars] \t\t" << comment << "\n";
  if (dg_parse == nullptr) {
    return;
  }

  DECL_IND2(ind)

  dump_sym_table(dg_parse->reg_sym_table, "reg_sym_table", "", ind2);

  //cout << ind2 << "root_object: \t\t" << dg_parse->root_object << " [object*] \t\t  object *root_instance;\n";
  dump_object_list(dg_parse->root_object, "root_object", "", ind2);

  //cout << ind2 << "current_object: \t\t" << dg_parse->current_object << " [object*] \t\t  region *current_region;\n";
  dump_object_list(dg_parse->current_object, "current_object", "", ind2);

  //cout << ind2 << "object_name_list: \t\t" << dg_parse->object_name_list << " [name_list*] \t\t  name_list *object_name_list_end;\n";
  dump_name_list(dg_parse->object_name_list, "object_name_list", "", ind2);

  cout << ind2 << "curr_file: \t\t" << str(dg_parse->curr_file) << " [char*] \t\t /* Name of MDL file currently being parsed */\n";

  dump_int_array(MAX_INCLUDE_DEPTH, "MAX_INCLUDE_DEPTH", (int*)dg_parse->line_num, "line_num", "/* Line numbers and filenames for all of the currently parsing files */", ind2);
  dump_string_array(MAX_INCLUDE_DEPTH, "MAX_INCLUDE_DEPTH", dg_parse->include_filename, "include_filename", "/* Line numbers and filenames for all of the currently parsing files */", ind2);

  cout << ind2 << "include_stack_ptr: \t\t" << dg_parse->include_stack_ptr << " [u_int] \t\t /* Stack pointer for filename/line number stack */\n";
  cout << ind2 << "comment_started: \t\t" << dg_parse->comment_started << " [int] \t\t /* Line number where last top-level (i.e. non-nested) multi-line (C-style) comment was started in the current MDL file. */\n";
}


void dump_dg_time_filename_list(dg_time_filename* fn, const char* name, const char* comment, const char* ind) {
  cout << ind << name << ": *\t\t" << (void*)fn << " [dyngeom_parse_vars] \t\t" << comment << "\n";
  dg_time_filename* curr = fn;
  int i = 0;
  while (curr != NULL) {
    cout << ind << i << ": " << curr->event_time << " - " << str(curr->mdl_file_path) << "\n";
    i++;
    curr = curr->next;
  }
}

/*extern "C"*/ void dump_volume(struct volume* s, const char* comment, unsigned int selected_details /* mask */) {

  cout << "********* volume dump :" << comment << "************ (START)\n";
  cout << comment << "\n";

  cout << "bond_angle: \t\t" << s->bond_angle << " [double] \t\t/* Defines the default bond rotation angle between molecules. Default is 0. */\n";

  // These are only used with dynamic geometry
  //cout << "dg_parse: *\t\t" << (void*)s->dg_parse << " [dyngeom_parse_vars] \n";
  dump_dyngeom_parse_vars(s->dg_parse, "dg_parse", "", "");

  cout << "dynamic_geometry_filename: *\t\t" << str(s->dynamic_geometry_filename) << " [char] \n";
  dump_molecules(s->num_all_molecules, s->all_molecules);

  cout << "names_to_ignore: *\t\t" << (void*)s->names_to_ignore << " [string_buffer] \n";

  /* Coarse partitions are input by the user */
  /* They may also be generated automagically */
  /* They mark the positions of initial partition boundaries */
  dump_double_array(s->nx_parts, "nx_parts", s->x_partitions, "x_partitions", "/* Coarse X partition boundaries */", "", 1000);
  dump_double_array(s->ny_parts, "ny_parts", s->y_partitions, "y_partitions", "/* Coarse Y partition boundaries */", "", 1000);
  dump_double_array(s->nz_parts, "nz_parts", s->z_partitions, "z_partitions", "/* Coarse Z partition boundaries */", "", 1000);

  cout << "mem_part_x: \t\t" << s->mem_part_x << " [int] \t\t/* Granularity of memory-partition binning for the X-axis */\n";
  cout << "mem_part_y: \t\t" << s->mem_part_y << " [int] \t\t/* Granularity of memory-partition binning for the Y-axis */\n";
  cout << "mem_part_z: \t\t" << s->mem_part_z << " [int] \t\t/* Granularity of memory-partition binning for the Z-axis */\n";
  cout << "mem_part_pool: \t\t" << s->mem_part_pool << " [int] \t\t/* Scaling factor for sizes of memory pools in each storage */\n";

  /* Fine partitions are intended to allow subdivision of coarse partitions */
  /* Subdivision is not yet implemented */
  dump_double_array(s->n_fineparts, "n_fineparts", s->x_fineparts, "x_fineparts", "/* Fine X partition boundaries */", "");
  dump_double_array(s->n_fineparts, "n_fineparts", s->y_fineparts, "y_fineparts", "/* Fine Y partition boundaries */", "");
  dump_double_array(s->n_fineparts, "n_fineparts", s->z_fineparts, "z_fineparts", "/* Fine Z partition boundaries */", "");

  cout << "periodic_traditional: \t\t" << s->periodic_traditional << " [bool] \n";

#ifdef DUMP_WAYPOINTS
  dump_waypoints(s->n_waypoints, s->waypoints);
#endif

  cout << "place_waypoints_flag: \t\t" << (int)s->place_waypoints_flag << " [byte] \t\t/* Used to save memory if waypoints not needed */\n";

#ifdef DUMP_SUBVOLUMES
  dump_subvolumes(s->n_subvols, "n_subvols", "/* How many coarse subvolumes? */", s->subvol, "subvol", "/* Array containing all subvolumes */", "");
#endif

  cout << "n_walls: \t\t" << s->n_walls << " [int] \t\t/* Total number of walls */\n";

  cout << "n_verts: \t\t" << s->n_verts << " [int] \t\t/* Total number of vertices */\n";

  dump_vector3_array(s->n_verts, "n_verts", s->all_vertices, "all_vertices", "/* Central repository of vertices with a partial order imposed by natural ordering of storages*/", "");

#ifdef DUMP_WALLS
  dump_wall_list_array(s->n_verts, "n_verts", s->walls_using_vertex, "walls_using_vertex", "/* Array of linked lists of walls using a vertex (has the size of all_vertices array) */", "");
#endif

  cout << "rx_hashsize: \t\t" << s->rx_hashsize << " [int] \t\t/* How many slots in our reaction hash table? */\n";
  cout << "n_reactions: \t\t" << s->n_reactions << " [int] \t\t/* How many reactions are there, total? */\n";


  cout << "reaction_hash: **\t\t" << (void**)s->reaction_hash << " [rxn] \t\t/* A hash table of all reactions. */\n";
  dump_reaction_hash_table(s->rx_hashsize, "rx_hashsize", s->reaction_hash, "reaction_hash", "/* A hash table of all reactions. */", "");

  cout << "tv_rxn_mem: *\t\t" << (void*)s->tv_rxn_mem << " [mem_helper] \t\t/* Memory to store time-varying reactions */\n";

  cout << "count_hashmask: \t\t" << s->count_hashmask << " [int] \t\t/* Mask for looking up count hash table */\n";
  cout << "count_hash: **\t\t" << (void**)s->count_hash << " [counter] \t\t/* Count hash table */\n";
  dump_schedule_helper(s->count_scheduler, "count_scheduler", "// When to generate reaction output", "", false);

  //cout << "counter_by_name: *\t\t" << (void*)s->counter_by_name << " [sym_table_head] \n";
  dump_sym_table(s->counter_by_name, "counter_by_name", "", "");

  dump_schedule_helper(s->volume_output_scheduler, "volume_output_scheduler", "/* When to generate volume output */", "", false);


  cout << "n_species: \t\t" << s->n_species << " [int] \t\t/* How many different species (molecules)? */\n";

  //NFSim stats
  cout << "n_NFSimSpecies: \t\t" << s->n_NFSimSpecies << " [int] \t\t/* number of graph patterns encountered during the NFSim simulation */\n";
  cout << "n_NFSimReactions: \t\t" << s->n_NFSimReactions << " [int] \t\t/* number of reaction rules discovered through an NFSim simulation */\n";
  cout << "n_NFSimPReactions: \t\t" << s->n_NFSimPReactions << " [int] \t\t/* number of potential reactions found in the reaction network (not necessarely triggered) */\n";

  cout << "species_list: **\t\t" << (void**)s->species_list << " [species] \t\t/* Array of all species (molecules). */\n";
  dump_species_list(s->n_species, "n_species", s->species_list, "species_list", "/* Array of all species (molecules). */", "");

  // dynamic geometry, not so important for now
  cout << "dynamic_geometry_flag: \t\t" << s->dynamic_geometry_flag << " [int] \t\t/* This is used to skip over certain sections in the parser when using dynamic geometries.*/ \n";
  cout << "disable_polygon_objects: \t\t" << s->disable_polygon_objects << " [int] \t\t/* This is used to skip over certain sections in the parser when using dynamic geometries.*/\n";


  //cout << "dynamic_geometry_head: *\t\t" << (void*)s->dynamic_geometry_head << " [dg_time_filename] \t\t/*List of all the dynamic geometry events that need to be scheduled*/\n";
  dump_dg_time_filename_list(s->dynamic_geometry_head, "dynamic_geometry_head", "/*List of all the dynamic geometry events that need to be scheduled*/", "  ");

  cout << "dynamic_geometry_events_mem: *\t\t" << (void*)s->dynamic_geometry_events_mem << " [mem_helper] \t\t/*Memory to store time and MDL names for dynamic geometry*/\n";

  //cout << "dynamic_geometry_scheduler: *\t\t" << (void*)s->dynamic_geometry_scheduler << " [schedule_helper] \t\t/*Scheduler for dynamic geometry (LATER)*/\n";
  dump_schedule_helper(s->dynamic_geometry_scheduler, "dynamic_geometry_scheduler", "/*Scheduler for dynamic geometry (LATER)*/", "", false);

  dump_schedule_helper(s->releaser, "releaser", "/* Scheduler for release events */", "", false);

  cout << "storage_allocator: *\t\t" << (void*)s->storage_allocator << " [mem_helper] \t\t/* Memory for storage list */\n";
  cout << "storage_head: *\t\t" << (void*)s->storage_head << " [storage_list] \t\t/* Linked list of all local  memory/schedulers\n";

  cout << "current_mol_id: \t\t" << s->current_mol_id << " [u_long] \t\t/* next unique molecule id to use*/\n";

  cout << "speed_limit: \t\t" << s->speed_limit << " [double] // How far can the fastest particle get in one timestep?\n";

  // symbol tables
  cout << "fstream_sym_table: *\t\t" << (void*)s->fstream_sym_table << " [sym_table_head] \t\t/* Global MDL file stream symbol hash table */\n";
  cout << "var_sym_table: *\t\t" << (void*)s->var_sym_table << " [sym_table_head] \t\t/* Global MDL variables symbol hash table */\n";
  cout << "rxn_sym_table: *\t\t" << (void*)s->rxn_sym_table << " [sym_table_head] \t\t/* RXN symbol hash table */\n";
  cout << "obj_sym_table: *\t\t" << (void*)s->obj_sym_table << " [sym_table_head] \t\t/* Objects symbol hash table */\n";
  cout << "reg_sym_table: *\t\t" << (void*)s->reg_sym_table << " [sym_table_head] \t\t/* Regions symbol hash table */\n";
  cout << "mol_sym_table: *\t\t" << (void*)s->mol_sym_table << " [sym_table_head] \t\t/* Molecule type symbol hash table */\n";
  cout << "rpat_sym_table: *\t\t" << (void*)s->rpat_sym_table << " [sym_table_head] \t\t/* Release pattern hash table */\n";
  cout << "rxpn_sym_table: *\t\t" << (void*)s->rxpn_sym_table << " [sym_table_head] \t\t/* Named reaction pathway hash table */\n";
  cout << "mol_ss_sym_table: *\t\t" << (void*)s->mol_ss_sym_table << " [sym_table_head] \t\t/* Spatially structured molecule symbol hash table */\n";

  // ???
  //cout << "root_object: *\t\t" << (void*)s->root_object << " [object] \t\t\n";
  dump_object_list(s->root_object, "root_object", "/* Root of the object template tree */", "");

  //cout << "root_instance: *\t\t" << (void*)s->root_instance << " [object] \t\t/* Root of the instantiated object tree */\n";
  dump_object_list(s->root_instance, "root_instance", "/* Root of the instantiated object tree */", "");

  cout << "periodic_box_obj: *\t\t" << (void*)s->periodic_box_obj << " [object] \n";


  // TODO3
  cout << "default_release_pattern: *\t\t" << (void*)s->default_release_pattern << " [release_pattern] \t\t/* release once at t=0 */\n";

  // ???
  cout << "volume_output_head: *\t\t" << (void*)s->volume_output_head << " [volume_output_item] \t\t/* List of all volume data output items */\n";


  //cout << "output_block_head: *\t\t" << (void*)s->output_block_head << " [output_block] \t\t/* Global list of reaction data output blocks */\n";
  dump_output_blocks(s->output_block_head, "output_block_head", "/* Global list of reaction data output blocks */", "");

  //cout << "output_request_head: *\t\t" << (void*)s->output_request_head << " [output_request] \t\t/* Global list linking COUNT statements to internal variables \n";
  dump_output_requests(s->output_request_head, "output_request_head", "/* Global list linking COUNT statements to internal variables*/", "");

  // some memories
  cout << "oexpr_mem: *\t\t" << (void*)s->oexpr_mem << " [mem_helper] \t\t/* Memory to store output_expressions */\n";
  cout << "outp_request_mem: *\t\t" << (void*)s->outp_request_mem << " [mem_helper] \t\t/* Memory to store output_requests */\n";
  cout << "counter_mem: *\t\t" << (void*)s->counter_mem << " [mem_helper] \t\t/* Memory to store counters (for counting molecules/reactions on regions)\n";
  cout << "trig_request_mem: *\t\t" << (void*)s->trig_request_mem << " [mem_helper] \t\t/* Memory to store listeners for trigger events */\n";
  cout << "magic_mem: *\t\t" << (void*)s->magic_mem << " [mem_helper] \t\t/* Memory used to store magic lists for reaction-triggered releases and such \n";

  //
  cout << "elapsed_time: \t\t" << s->elapsed_time << " [double] \t\t/* Used for concentration measurement */\n";

  /* Visualization state */
  dump_viz_output_block(s->viz_blocks, "viz_blocks", "/* VIZ_OUTPUT blocks from file */", "");

  dump_species(s->all_mols, "all_mols", "/* Refers to ALL_MOLECULES keyword */", "");
  dump_species(s->all_volume_mols, "all_volume_mols", "/* Refers to ALL_VOLUME_MOLECULES keyword */", "");
  dump_species(s->all_surface_mols, "all_surface_mols", "/* Refers to ALL_SURFACE_MOLECULES keyword */", "");

  cout << "time_unit: \t\t" << s->time_unit << " [double] \t\t/* Duration of one global time step in real time, Used to convert between real time and internal time */\n";
  cout << "time_step_max: \t\t" << s->time_step_max << " [double] \t\t/* Maximum internal time that a molecule may diffuse */\n";

  cout << "grid_density: \t\t" << s->grid_density << " [[double] \t\t/* Density of grid for surface molecules, number per um^2 */\n";
  cout << "length_unit: \t\t" << s->length_unit << " [double] \t\t/* Internal unit of distance, 1/sqrt(grid_density), in microns\n";

  cout << "r_length_unit: \t\t" << s->r_length_unit << " [double] \t\t/* Reciprocal of length_unit to avoid division */\n";
  cout << "rx_radius_3d: \t\t" << s->rx_radius_3d << " [double] \t\t/* Interaction radius for reactions between volume molecules \n";

  cout << "space_step: \t\t" << s->space_step << " [double] \t\t/* User-supplied desired average diffusion distance for volume molecules\n";

  cout << "r_step: *\t\t" << (void*)s->r_step << " [double] \t\t/* Lookup table of 3D diffusion step lengths */\n";
  cout << "d_step: *\t\t" << (void*)s->d_step << " [double] \t\t/* Lookup table of 3D diffusion direction vectors */\n";
  cout << "r_step_surface: *\t\t" << (void*)s->r_step_surface << " [double] \t\t/* Lookup table of 2D diffusion step lengths */\n";
  cout << "r_step_release: *\t\t" << (void*)s->r_step_release << " [double] \t\t/* Lookup table of diffusion lengths for 3D release */\n";
  cout << "radial_subdivisions: \t\t" << s->radial_subdivisions << " [u_int] \t\t/* Size of 2D and 3D step length lookup tables */\n";
  cout << "radial_directions: \t\t" << s->radial_directions << " [u_int] \t\t/* Requested size of 3D direction lookup table */\n";
  cout << "num_directions: \t\t" << s->num_directions << " [u_int] \t\t/* Actual size of 3D direction lookup table */\n";
  cout << "directions_mask: \t\t" << s->directions_mask << " [int] \t\t/* Mask to obtain RNG bits for direction lookup */\n";
  cout << "fully_random: \t\t" << s->fully_random << " [int] \t\t/* If set, generate directions with trig functions instead of lookup table \n";

  cout << "dissociation_index: \t\t" << s->dissociation_index << " [int] \t\t/* Used to keep 3D products from reacting with each  other too soon \n";

  // chkpt
  cout << "chkpt_iterations: \t\t" << s->chkpt_iterations << " [long long] \t\t/* Number of iterations to advance before checkpointing\n";
  cout << "chkpt_init: \t\t" << s->chkpt_init << " [u_int] \t\t/* Set if this is the initial run of a simulation with no previous checkpoints\n";
  cout << "chkpt_flag: \t\t" << s->chkpt_flag << " [u_int] \t\t/* Set if there are any CHECKPOINT statements in mdl file */\n";
  cout << "chkpt_seq_num: \t\t" << s->chkpt_seq_num << " [u_int] \t\t/* Number of current run in checkpoint sequence */\n";
  cout << "keep_chkpts: \t\t" << s->keep_chkpts << " [int] \t\t/* flag to indicate if checkpoints should be kept */\n";
  cout << "chkpt_infile: *\t\t" << (void*)s->chkpt_infile << " [char] \t\t/* Name of checkpoint file to read from */\n";
  cout << "chkpt_outfile: *\t\t" << (void*)s->chkpt_outfile << " [char] \t\t/* Name of checkpoint file to write to */\n";
  cout << "chkpt_byte_order_mismatch: \t\t" << s->chkpt_byte_order_mismatch << " [u_int] \t\t/* Flag that defines whether mismatch in byte order exists between the saved checkpoint file and the machine reading it\n";
  cout << "chkpt_start_time_seconds: \t\t" << s->chkpt_start_time_seconds << " [double] \t\t/* start of the simulation time (in sec) for new checkpoint\n";

  cout << "current_time_seconds: \t\t" << s->current_time_seconds << " [double] \t\t/* current simulation time in seconds */\n";
  /* simulation start time (in seconds) or time of most recent checkpoint */
  cout << "simulation_start_seconds: \t\t" << s->simulation_start_seconds << " [double] \n";

  cout << "diffusion_number: \t\t" << s->diffusion_number << " [long long] \t\t/* Total number of times molecules have had their positions updated */\n";

  cout << "diffusion_cumtime: \t\t" << s->diffusion_cumtime << " [double] \t\t/* Total time spent diffusing by all molecules */\n";
  cout << "ray_voxel_tests: \t\t" << s->ray_voxel_tests << " [long long] \t\t/* How many ray-subvolume intersection tests have we performed */\n";

  cout << "ray_polygon_tests: \t\t" << s->ray_polygon_tests << " [long long] \t\t/* How many ray-polygon intersection tests have  we performed */\n";

  cout << "ray_polygon_colls: \t\t" << s->ray_polygon_colls << " [long long] \t\t/* How many ray-polygon intersections have occured */\n";

  cout << "dyngeom_molec_displacements: \t\t" << s->dyngeom_molec_displacements << " [long long] \t\t/* Total number of dynamic geometry molecule displacements */\n";

  /* below "vol" means volume molecule, "surf" means surface molecule */
  cout << "vol_vol_colls: \t\t" << s->vol_vol_colls << " [long long] \t\t/* How many vol-vol collisions have occured */\n";
  cout << "vol_surf_colls: \t\t" << s->vol_surf_colls << " [long long] \t\t/* How many vol-surf collisions have occured */\n";
  cout << "surf_surf_colls: \t\t" << s->surf_surf_colls << " [long long] \t\t/* How many surf-surf collisions have occured */\n";
  cout << "vol_wall_colls: \t\t" << s->vol_wall_colls << " [long long] \t\t/* How many vol-wall collisions have occured */\n";
  cout << "vol_vol_vol_colls: \t\t" << s->vol_vol_vol_colls << " [long long] // How many vol-vol-vol collisions have occured\n";
  cout << "vol_vol_surf_colls: \t\t" << s->vol_vol_surf_colls << " [long long] \t\t/* How many vol-vol-surf collisions have occured */\n";
  cout << "vol_surf_surf_colls: \t\t" << s->vol_surf_surf_colls << " [long long] \t\t/* How many vol-surf-surf collisions have occured */\n";

  cout << "surf_surf_surf_colls: \t\t" << s->surf_surf_surf_colls << " [long long] \t\t/* How many surf-surf-surf collisions have occured */\n";


  cout << "bb_llf: \t\t" << s->bb_llf << " [vector3] \t\t/* llf corner of world bounding box */\n";
  cout << "bb_urb: \t\t" << s->bb_urb << " [vector3] \t\t/* urb corner of world bounding box */\n";

  cout << "rng: *\t\t" << (void*)s->rng << " [rng_state] \t\t/* State of the random number generator (currently isaac64) */\n";

  cout << "init_seed: \t\t" << s->init_seed << " [u_int] \t\t/* Initial seed value for random number generator */\n";

  cout << "current_iterations: \t\t" << s->current_iterations << " [long long] \t\t/* How many iterations have been run so far */\n";
  // Starting iteration number for current run or iteration of most recent
  // checkpoint
  cout << "start_iterations: \t\t" << s->start_iterations << " [long long] \n";
  cout << "last_timing_time: \t\t" << s->last_timing_time << " [timeval] \t\t/* time and iteration of last timing event */\n";
  cout << "last_timing_iteration: \t\t" << s->last_timing_iteration << " [long long] \t\t/* during the main run_iteration loop */\n";

  cout << "procnum: \t\t" << s->procnum << " [int] \t\t/* Processor number for a parallel run */\n";
  cout << "quiet_flag: \t\t" << s->quiet_flag << " [int] \t\t/* Quiet mode */\n";
  cout << "with_checks_flag: \t\t" << s->with_checks_flag << " [int] \t\t/* Check geometry for overlapped walls? */\n";

  cout << "coll_mem: *\t\t" << (void*)s->coll_mem << " [mem_helper] \t\t/* Collision list */\n";
  cout << "sp_coll_mem: *\t\t" << (void*)s->sp_coll_mem << " [mem_helper] \t\t/* Collision list (trimol) */\n";
  cout << "tri_coll_mem: *\t\t" << (void*)s->tri_coll_mem << " [mem_helper] \t\t/* Collision list (trimol) */\n";
  cout << "exdv_mem: *\t\t" << (void*)s->exdv_mem << " [mem_helper] // Vertex lists f    ct interaction disk area\n";

  /* Current version number. Format is "3.XX.YY" where XX is major release
   * number (for new features) and YY is minor release number (for patches) */
  cout << "mcell_version: *\t\t" << (void*)s->mcell_version << " [const char] \n";

  cout << "use_expanded_list: \t\t" << s->use_expanded_list << " [int] \t\t/* If set, check neighboring subvolumes for mol-mol interactions */\n";

  cout << "randomize_smol_pos: \t\t" << s->randomize_smol_pos << " [int] \t\t/* If set, always place surface molecule at random location instead of center of grid */\n";

  cout << "vacancy_search_dist2: \t\t" << s->vacancy_search_dist2 << " [double] \t\t/* Square of distance to search for free grid  location to place surface product */\n";

  cout << "surface_reversibility: \t\t" << s->surface_reversibility << " [byte] \t\t/* If set, match unbinding diffusion distribution  to binding distribution at surface */\n";

  cout << "volume_reversibility: \t\t" << s->volume_reversibility << " [byte] \t\t/* If set, match unbinding diffusion distribution to binding distribution in volume */\n";


  /* If set to NEAREST_TRIANGLE, molecules are moved to a random location
   * slightly offset from the enclosing wall. If set to NEAREST_POINT, then
   * they are moved to the closest point on that wall (still slightly offset).
   * */
  cout << "dynamic_geometry_molecule_placement: \t\t" << s->dynamic_geometry_molecule_placement << " [int] \n";

  /* MCell startup command line arguments */
  cout << "seed_seq: \t\t" << s->seed_seq << " [u_int] \t\t/* Seed for random number generator */\n";
  cout << "iterations: \t\t" << s->iterations << " [long long] \t\t/* How many iterations to run */\n";
  cout << "log_freq: \t\t" << s->log_freq << " [unsigned long] \t\t/* Interval between simulation progress reports, default scales as sqrt(iterations) */\n";

  cout << "mdl_infile_name: *\t\t" << (void*)s->mdl_infile_name << " [char] \t\t/* Name of MDL file specified on command line */\n";
  cout << "curr_file: *\t\t" << (void*)s->curr_file << " [const char ] \t\t/* Name of MDL file currently being parsed */\n";

  // XXX: Why do we allocate this on the heap rather than including it inline?
  cout << "notify: *\t\t" << (void*)s->notify << " [notifications] \t\t/* Notification/warning/output flags */\n";

  cout << "clamp_list: *\t\t" << (void*)s->clamp_list << " [ccn_clamp_data] \t\t/* List of objects at which volume molecule concentrations should be clamped */\n";



  /* Flags for asynchronously-triggered checkpoints */

  /* Flag indicating whether a checkpoint has been requested. */
  cout << "checkpoint_requested: \t\t" << s->checkpoint_requested << " [checkpoint_request_type_t] \n";
  cout << "checkpoint_alarm_time: \t\t" << s->checkpoint_alarm_time << " [unsigned int] \t\t// number of seconds between checkpoints\n";
  cout << "continue_after_checkpoint: \t\t" << s->continue_after_checkpoint << " [int] \t\t/* 0: exit after chkpt, 1: continue after chkpt */\n";
  cout << "last_checkpoint_iteration: \t\t" << s->last_checkpoint_iteration << " [long long] \t\t/* Last iteration when chkpt was created */\n";
  cout << "begin_timestamp: \t\t" << s->begin_timestamp << " [time_t] \t\t/* Time since epoch at beginning of 'main' */\n";
  cout << "initialization_state: *\t\t" << (void*)s->initialization_state << " [char] \t\t/* NULL after initialization completes */\n";
  cout << "rxn_flags: \t\t" << s->rxn_flags << " [reaction_flags] \n";
  /* shared walls information per mesh vertex is created when there are
   reactions present with more than one surface reactant or more than one
   surface product */
  cout << "create_shared_walls_info_flag: \t\t" << s->create_shared_walls_info_flag << " [int] \n";
  /* resource usage during initialization */
  cout << "u_init_time: \t\t" << s->u_init_time << " [timeval] \t\t/* user time */\n";
  cout << "s_init_time: \t\t" << s->s_init_time << " [timeval] \t\t/* system time */\n";
  cout << "t_start: \t\t" << s->t_start << " [time_t] \t\t/* global start time */\n";
  cout << "reaction_prob_limit_flag: \t\t" << s->reaction_prob_limit_flag << " [byte] \t\t/* checks whether there is at least one reaction with probability greater than 1 including variable rate reactions */\n";



  //JJT: Checks if we will be communicating with nfsim
  cout << "nfsim_flag: \t\t" << s->nfsim_flag << " [int] \n";
  cout << "global_nfsim_volume: *\t\t" << (void*)s->global_nfsim_volume << " [species] \n";
  cout << "global_nfsim_surface: *\t\t" << (void*)s->global_nfsim_surface << " [species] \n";

  cout << "species_mesh_transp: *\t\t" << (void*)s->species_mesh_transp << " [pointer_hash] \n";

  cout << "********* volume dump :" << comment << "************ (END)\n";

  cout.flush();
}


// ----------- other debug functions -------

string collision_flags_to_str(int flags) {
  string res;
  if (flags == COLLIDE_MISS) {
    res = "COLLIDE_MISS";
    return res;
  }
  if (flags == COLLIDE_REDO) {
    res = "COLLIDE_REDO";
    return res;
  }

#define CHECK_FLAG(F) if (flags & F) res += #F ","

  CHECK_FLAG(COLLIDE_FRONT);
  CHECK_FLAG(COLLIDE_BACK);
  CHECK_FLAG(COLLIDE_VOL_M);
  CHECK_FLAG(COLLIDE_SV_NX);
  CHECK_FLAG(COLLIDE_SV_PX);
  CHECK_FLAG(COLLIDE_SV_NY);
  CHECK_FLAG(COLLIDE_SV_PY);
  CHECK_FLAG(COLLIDE_SV_NZ);
  CHECK_FLAG(COLLIDE_SV_PZ);
  CHECK_FLAG(COLLIDE_MASK);
  CHECK_FLAG(COLLIDE_WALL);
  CHECK_FLAG(COLLIDE_VOL);
  CHECK_FLAG(COLLIDE_SUBVOL);
  CHECK_FLAG(COLLIDE_VOL_VOL);
  CHECK_FLAG(COLLIDE_VOL_SURF);
  CHECK_FLAG(COLLIDE_SURF_SURF);
  CHECK_FLAG(COLLIDE_SURF);

#undef CHECK_FLAG
  return res;
}

void dump_one_collision(collision* col) {
  cout << "next: *\t\t" << col->next << " [collision] \t\t\n";
  cout << "t: \t\t" << col->t << " [double] \t\t/* Time of collision (may be slightly early) */\n";

  cout << "target: *\t\t" << (void*)col->target << " [void] \t\t/* Thing that we hit: wall, molecule, subvol etc */\n";
  cout << "what: \t\t" << collision_flags_to_str(col->what) << " [int] \t\t/* Target-type Flags: what kind of thing did we hit? */\n";
  cout << "intermediate: *\t\t" << col->intermediate << " [rxn] \t\t/* Reaction that told us we could hit it */\n";
  //if (col->intermediate != nullptr && col->intermediate->sym != nullptr && col->intermediate->sym->name != nullptr) {
  //  cout << "intermediate->sym->name: *\t\t" << col->intermediate->sym->name << "\n";
  //}
  cout << "loc: \t\t" << col->loc << " [vector3] \t\t/* Location of impact */\n";
}

void dump_collisions(collision* shead) {
  /*cout << "Collision list: " << ((shead == nullptr) ? "EMPTY" : "") << "\n";
*/
  int i = 0;
  collision* ptr = shead;
  while (ptr != NULL) {
    if (ptr->what & COLLIDE_VOL /* != 0 && ptr->t < 1.0 && ptr->t >= 0.0*/) {
      cout << "  " << "mol collision " << i << ": "
          //<< "diff_idx: " << ptr-> diffused_molecule_idx
          << "coll_idx: " << ((volume_molecule*)ptr->target)->id
          << ", time: " << ptr->t
          << ", pos: " << ptr->loc
          << "\n";
      i++;
    }
    else if (
        (ptr->what & COLLIDE_SUBVOL) == 0
        && (ptr->what == COLLIDE_REDO
            || (ptr->what & COLLIDE_FRONT) != 0
            || (ptr->what & COLLIDE_BACK) != 0
           )
    ) {
#ifndef NODEBUG_WALL_COLLISIONS
      const char* name = ((wall*)ptr->target)->parent_object->sym->name;
      cout << "  " << "wall collision " << i << ": "
          //<< "diff_idx: " << ptr-> diffused_molecule_idx
          << "obj name: " << ((name != nullptr) ? name : "")
          << ", wall side: " << ((wall*)ptr->target)->side
          << ", time: " << ptr->t
          << ", pos: " << ptr->loc
          << "\n";
#endif
      i++;
    }
    ptr = ptr->next;
  }
}

// -----------  differential dumps ---------------

string get_species_name(volume_molecule* vm) {
  return vm->properties->sym->name;
}

string get_species_name(surface_molecule* sm) {
  return sm->properties->sym->name;
}

void dump_volume_molecule(
    volume_molecule* vm,
    const char* ind,
    bool for_diff,
    const char* extra_comment,
    unsigned long long iteration,
    double time,
    bool print_position
) {
  if (!for_diff) {
    cout << ind << "id: \t\t" << vm->id << " [u_long] \t\t\n";
    cout << ind << "pos: \t\t" << vm->pos << " [vector3] \t\t/* Position in space */\n";
    cout << ind << "properties: *\t\t" << vm->properties << " [species] \t\t\n";
    cout << ind << "  species name: *\t\t" << vm->properties->sym->name << " [char] \t\t\n";
  }
  else {
    cout << extra_comment << "it:" << iteration << ", idx:" << vm->id;

    if ((vm->properties->flags & EXTERNAL_SPECIES) == 0) {
      cout << ", species: " << get_species_name(vm);
    }
    else {
      cout << ", species: " << graph_pattern_to_bngl(vm->graph_data->graph_pattern);
    }
    if (print_position) {
      cout << ", pos:" << vm->pos;
    }
    cout << ", flags:" << get_molecule_flags_string(vm->flags, false) << ", time: " << time << "\n";
  }
}

void dump_surface_molecule(
    surface_molecule* sm,
    const char* ind,
    bool for_diff,
    const char* extra_comment,
    unsigned long long iteration,
    double time,
    bool print_position
) {
  if (!for_diff) {
    cout << ind << "id: \t\t" << sm->id << " [u_long] \t\t\n";
    cout << ind << "pos: \t\t" << sm->s_pos << " [vector2] \t\t/* Position in space */\n";
    cout << ind << "properties: *\t\t" << sm->properties << " [species] \t\t\n";
    cout << ind << "  species name: *\t\t" << sm->properties->sym->name << " [char] \t\t\n";
  }
  else {
    cout << extra_comment << "it:" << iteration << ", idx:" << sm->id
        << ", species: " << get_species_name(sm);
    if (print_position) {
        cout << ", pos:" << sm->s_pos
            << ", orient:" << sm->orient
            << ", wall side: " << ((sm->grid == nullptr || sm->grid->surface == nullptr) ? -1 : sm->grid->surface->side)
            << ", grid index: " << sm->grid_index;
    }
    cout << ", flags:" << get_molecule_flags_string(sm->flags, false) << ", time: " << time << "\n";
  }
}

void dump_vector2(struct vector2 vec, const char* extra_comment) {
  cout << extra_comment << vec << "\n";
}

void dump_vector3(struct vector3 vec, const char* extra_comment) {
  cout << extra_comment << vec << "\n";
}

void dump_tile_neighbors_list(struct tile_neighbor *tile_nbr_head, const char* extra_comment, const char* ind) {
  struct tile_neighbor *curr = tile_nbr_head;
  int i = 0;
  std::cout << ind << extra_comment << "\n";
  while (curr != nullptr) {
    std::cout << ind << i << ": " << curr->idx << "\n";
    curr = curr->next;
    i++;
  }
}

void dump_processing_reaction(
    long long it,
    struct vector3 *hitpt, double t,
    struct rxn *rx, /*int path,*/
    struct abstract_molecule *reacA,
    struct abstract_molecule *reacB,
    struct wall *w
) {
  assert(reacA != nullptr);
  bool two_reactants = rx->n_reactants == 2;

  cout << "Processing reaction:it:" << it << ", ";

  if (two_reactants) {

    if (reacB != nullptr) {
      cout <<
        "bimol rxn" <<
        ", idA:"  << reacA->id <<
        ", idB:"  << reacB->id <<
        //TODO ", rxn: " << rx->to_string(p) <<
        ", time: " << t;
      if (hitpt != nullptr) {
        cout << ", pos " << *hitpt;
      }
    }
    else {
      assert(w != nullptr);
      cout <<
        "wall collision" <<
        ", idA:"  << reacA->id <<
        //", wall_index:"  << w->side <<
        //TODO ", rxn: " << rx->to_string(p) <<
        ", time: " << t;

      if (hitpt != nullptr) {
        cout << ", pos " << *hitpt;
      }
    }
  }
  else {
    cout <<
      "unimol rxn" <<
      ", idA:"  << reacA->id <<
      // ", rxn: " << rx->to_string(p) <<
      ", time: " << t;
  }
  cout << "\n";
}
