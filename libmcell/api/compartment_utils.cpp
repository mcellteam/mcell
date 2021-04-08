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

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkSelectEnclosedPoints.h>

#include "compartment_utils.h"

#include "api/api_common.h"
#include "counted_volume_utils.h"
#include "geometry_object.h"

using namespace std;

namespace MCell {

using namespace CountedVolumeUtils;

namespace API {

// note: similar code is in convert_and_set_geometry_object_polydata but uses different input data
static void convert_compartment_objects_to_geom_object_infos(
    const std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects,
    GeomObjectInfoVector& compartment_infos) {
  bool res = true;

  for (const auto& obj_ptr: compartment_objects) {
    const API::GeometryObject& obj = *obj_ptr;

    // we need to convert each geometry object into VTK's polydata representation
    // example: https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/TriangleArea
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

    GeomObjectInfo obj_info(obj.name);
    obj_info.polydata = vtkSmartPointer<vtkPolyData>::New();

    // first collect vertices
    uint_set<vertex_index_t> vertex_indices;
    for (const std::vector<float_t>& v: obj.vertex_list) {
      assert(v.size() == 3);
      points->InsertNextPoint(v[0], v[1], v[2]);
    }
    // store triangles
    for (const std::vector<int>& wi: obj.wall_list) {
      assert(wi.size() == VERTICES_IN_TRIANGLE);

      vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

      for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
        triangle->GetPointIds()->SetId(i, wi[i]);
      }

      triangles->InsertNextCell(triangle);
    }

    // create input polydata
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetPolys(triangles);

    // clean them up
    vtkSmartPointer<vtkTriangleFilter> tri = vtkSmartPointer<vtkTriangleFilter>::New();
    tri->SetInputData(polydata);
    vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputConnection(tri->GetOutputPort());
    clean->Update();

    int closed = vtkSelectEnclosedPoints::IsSurfaceClosed(clean->GetOutput());
    if (closed != 1) {
      throw RuntimeError("Compartment object must be closed, error for " + obj.name + ".");
    }

    // and finally store the points and faces
    obj_info.polydata = clean->GetOutput();

    compartment_infos.push_back(obj_info);
  }
}


static std::shared_ptr<API::GeometryObject> get_obj_by_name(
    std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects,
    const string& name) {
  for (auto& obj: compartment_objects) {
    if (obj->name == name) {
      return obj;
    }
  }
  release_assert(false);
  return std::shared_ptr<API::GeometryObject>(nullptr);
}


void set_parent_and_children_compartments(
    std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects) {

#if 0
  auto obj_cp = std::shared_ptr<API::GeometryObject>(nullptr);
  for (auto obj: compartment_objects) {
    release_assert(obj->name == "CP" || obj->name == "EC"); // limited for now
    if (obj->name == "CP") {
      obj_cp = obj;
    }
  }

  auto obj_ec = std::shared_ptr<API::GeometryObject>(nullptr);
  for (auto obj: compartment_objects) {
    if (obj->name == "EC") {
      obj->child_compartments.insert(obj_cp);
      obj_ec = obj;
    }
  }
  if (is_set(obj_cp) && is_set(obj_ec)) {
    obj_cp->parent_compartment = obj_ec;
  }

  return true;
#endif

  GeomObjectInfoVector compartment_infos;
  convert_compartment_objects_to_geom_object_infos(compartment_objects, compartment_infos);

  // get information on what objects are contained within each other
  // e.g. CP -> {EC} - CP is contained in EC, mapping contains all compartments, not just the direct parent
  ContainmentMap contained_in_mapping;
  IntersectingSet intersecting_object_infos;
  bool ok = compute_containement_mapping(
      nullptr, compartment_infos, contained_in_mapping, intersecting_object_infos);
  if (!ok) {
    throw RuntimeError("Unexpected error while computing compartment hierarchy.");
  }

  if (!intersecting_object_infos.empty()) {
    string names;
    for (const GeomObjectInfo& info: intersecting_object_infos) {
      names += info.name + ", ";
    }
    names = names.substr(0, names.size() - 2); // remove last comma
    throw RuntimeError("Error while computing compartment hierarchy, following compartments intersect " +
        names + ".");
  }

  // set children and parents
  for (const auto& pair_name_parents: contained_in_mapping) {
    auto child = get_obj_by_name(compartment_objects, pair_name_parents.first.name);

    // select the direct parent
    const GeomObjectInfo* direct_parent_info =
        get_direct_parent_info(pair_name_parents.first, contained_in_mapping);

    if (direct_parent_info == nullptr) {
      // no parent, children set themselves as children of this compartment
    }
    else {
      auto parent = get_obj_by_name(compartment_objects, direct_parent_info->name);

      parent->child_compartments.insert(child);
      if (is_set(child->parent_compartment)) {
        // this should not really occur but lets better report it
        throw RuntimeError("Compartment object " + child->name + " cannot have multiple parent compartments " +
            child->parent_compartment->name + " and " + parent->name + ".");
      }
      child->parent_compartment = parent;
    }
  }
}


} // namespace API
} // namespace MCell
