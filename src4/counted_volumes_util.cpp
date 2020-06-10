/******************************************************************************
 *
 * Copyright (C) 2020 by
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

#include <vector>
#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkCollisionDetectionFilter.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkTransform.h>
#include <vtkPointData.h>

#include "logging.h"

#include "defines.h"
#include "world.h"
#include "geometry.h"

using namespace std;

namespace MCell {

namespace CountedVolumesUtil {

// FIXME: use directly GeometryObject where it will have the partition ID
struct GeomObjectInfo {

  GeomObjectInfo(const partition_id_t partition_id_, const geometry_object_id_t geometry_object_id_)
    : partition_id(partition_id_), geometry_object_id(geometry_object_id_) {
  }

  // we are using IDs because in the future we might need to create new
  // geometry objects and pointers would become invalidated
  partition_id_t partition_id;
  geometry_object_id_t geometry_object_id;

  vtkSmartPointer<vtkPolyData> polydata;


  GeometryObject& get_geometry_object_noconst(World* world) const {
    Partition& p = world->get_partition(partition_id);
    return p.get_geometry_object(geometry_object_id);
  }

  const GeometryObject& get_geometry_object(const World* world) const {
    const Partition& p = world->get_partition(partition_id);
    return p.get_geometry_object(geometry_object_id);
  }

  // comparison uses geometry_object_id only, used in maps
  bool operator < (const GeomObjectInfo& second) const {
    return geometry_object_id < second.geometry_object_id;
  }

  bool operator == (const GeomObjectInfo& second) const {
    return geometry_object_id == second.geometry_object_id;
  }

};


typedef vector<GeomObjectInfo> GeomObjectInfoVector;

// containment mapping of counted geometry objects
typedef std::map<GeomObjectInfo, uint_set<GeomObjectInfo>> ContainmentMap;
typedef uint_set<geometry_object_id_t> IntersectingSet;

enum class ContainmentResult {
  Error,
  Disjoint,
  Identical,
  Intersect,
  Obj1InObj2,
  Obj2InObj1
};



static bool convert_objects_to_clean_polydata(World* world, GeomObjectInfoVector& counted_objects) {
  bool res = true;

  for (GeomObjectInfo& obj_info: counted_objects) {

    Partition& p = world->get_partition(obj_info.partition_id);
    GeometryObject& obj = p.get_geometry_object(obj_info.geometry_object_id);

    // we need to convert each geometry object into VTK's polydata representation
    // example: https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/TriangleArea
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    obj_info.polydata = vtkSmartPointer<vtkPolyData>::New();

    // first collect vertices
    uint_set<vertex_index_t> vertex_indices;
    for (wall_index_t wi: obj.wall_indices) {
      const Wall& w = p.get_wall(wi);

      for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
        vertex_indices.insert(w.vertex_indices[i]);
      }
    }

    // store vertices and create mapping
    // mcell vertex index -> vtk vertex index
    map<vertex_index_t, uint> vertex_mapping;
    uint curr_vtk_index = 0;
    for (vertex_index_t vi: vertex_indices) {

      Vec3 pt = p.get_geometry_vertex(vi);
      points->InsertNextPoint(pt.x, pt.y, pt.z);

      vertex_mapping[vi] = curr_vtk_index;
      curr_vtk_index++;
    }

    // store triangles
    for (wall_index_t wi: obj.wall_indices) {
      const Wall& w = p.get_wall(wi);

      vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

      for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
        triangle->GetPointIds()->SetId(i, vertex_mapping[w.vertex_indices[i]]);
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
      mcell_warn("Counting object must be closed, error for %s.", obj.name.c_str());
      res = false;
      continue;
    }

    // and finally store the points and faces
    obj_info.polydata = clean->GetOutput();

    // also copy this information to our object
    //obj.counted_volume_polydata.TakeReference(clean->GetOutput());
    obj.counted_volume_polydata = clean->GetOutput();
  }

  return res;
}


// the objects do not collide
static bool is_noncolliding_obj1_fully_contained_in_obj2(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2) {
  // NOTE: it will be sufficient to check just one point whether it is outside since we know that there is no intersect,
  // optimize in the future

  // is the object contained?
  auto select_enclosed_points = vtkSmartPointer<vtkSelectEnclosedPoints>::New();

  select_enclosed_points->SetSurfaceData(poly2); // poly 2 is supposed to be the larger object
  select_enclosed_points->SetInputData(poly1); // and we are trying whether poly1 fits into that
  select_enclosed_points->Update();

  vtkDataArray* inside_array = vtkDataArray::SafeDownCast(
      select_enclosed_points->GetOutput()->GetPointData()->GetAbstractArray("SelectedPoints"));

  for(vtkIdType i = 0; i < inside_array->GetNumberOfTuples(); i++)
  {
    if(inside_array->GetComponent(i,0) == 1)
    {
      // there was no collision so if any of the points is inside, the object is inside
      return true;
    }
  }
  return false;
}


static bool objs_have_identical_points(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2) {
  vtkSmartPointer<vtkPoints> points1 = poly1->GetPoints();
  vtkSmartPointer<vtkPoints> points2 = poly2->GetPoints();

  vtkIdType num_points1 = points1->GetNumberOfPoints();
  vtkIdType num_points2 = points2->GetNumberOfPoints();

  if (num_points1 != num_points2) {
    return false;
  }

  bool all_same = true;

  double verts1[3];
  double verts2[3];

  for(vtkIdType i = 0; i < num_points1; i++)
  {
    points1->GetPoint(i, verts1);
    points2->GetPoint(i, verts2);

    if (verts1[0] != verts2[0] || verts1[1] != verts2[1] || verts1[2] != verts2[2])
    {
      return false;
    }
  }

  return all_same;
}


static ContainmentResult geom_object_containment_test(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2) {

  // counting objects must be closed (already checked during conversion)
  assert(vtkSelectEnclosedPoints::IsSurfaceClosed(poly1) == 1);
  assert(vtkSelectEnclosedPoints::IsSurfaceClosed(poly2) == 1);

  // 1) do they collide?
  auto matrix1 = vtkSmartPointer<vtkMatrix4x4>::New();
  auto transform0 = vtkSmartPointer<vtkTransform>::New();
  vtkSmartPointer<vtkCollisionDetectionFilter> collide = vtkSmartPointer<vtkCollisionDetectionFilter>::New();

  collide->SetInputData( 0, poly1);
  collide->SetTransform(0, transform0);

  collide->SetInputData( 1, poly2);
  collide->SetMatrix(1, matrix1);

  collide->SetBoxTolerance(0.0);
  collide->SetCellTolerance(0.0);
  collide->SetNumberOfCellsPerNode(2);

  collide->SetCollisionModeToFirstContact();

  collide->GenerateScalarsOn();
  collide->Update();

  if (collide->GetNumberOfContacts() == 0) {
    // objects do not collide

    if (is_noncolliding_obj1_fully_contained_in_obj2(poly1, poly2)) {
      return ContainmentResult::Obj1InObj2;
    }
    else if (is_noncolliding_obj1_fully_contained_in_obj2(poly2, poly1)) {
      return ContainmentResult::Obj2InObj1;
    }
    else {
      return ContainmentResult::Disjoint;
    }
  }
  else {
    // the objects collide

    if (objs_have_identical_points(poly1, poly2)) {
      // do not necessarily have to be identical, but let's assume that if the points are the same,
      // the objects are the same
      return ContainmentResult::Identical;
    }
    else {

      return ContainmentResult::Intersect;

      #if 0
        // this is how intersect is computed once it will be needed
        vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
          vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
        booleanOperation->SetOperationToIntersection();
        booleanOperation->SetInputData( 0, poly1 );
        booleanOperation->SetInputData( 1, poly2 );
        booleanOperation->Update();
      #endif
    }
  }
}


// returns false if theee was any error
static bool compute_containement_mapping(
    const World* world, const GeomObjectInfoVector& counted_objects,
    ContainmentMap& contained_in_mapping,
    IntersectingSet& intersecting_objects
) {

  // keep it simple for now, let's just compute 'contained in' relation for each pair of objects
  // can be optimized in the future

  bool res = true;

  contained_in_mapping.clear();
  for (uint obj1 = 0; obj1 < counted_objects.size(); obj1++) {
    for (uint obj2 = obj1 + 1; obj2 < counted_objects.size(); obj2++) {
      ContainmentResult containment_res =
          geom_object_containment_test(counted_objects[obj1].polydata, counted_objects[obj2].polydata);

      switch (containment_res) {
        case ContainmentResult::Obj1InObj2:
          contained_in_mapping[counted_objects[obj1]].insert(counted_objects[obj2]);
          break;

        case ContainmentResult::Obj2InObj1:
          contained_in_mapping[counted_objects[obj2]].insert(counted_objects[obj1]);
          break;

        case ContainmentResult::Disjoint:
          // nothing to do
          break;

        case ContainmentResult::Intersect:
          // we are not processing all pairs so we need to insert both object ids
          intersecting_objects.insert(counted_objects[obj1].geometry_object_id);
          intersecting_objects.insert(counted_objects[obj2].geometry_object_id);
          break;

        case ContainmentResult::Identical:
        case ContainmentResult::Error: {
            std::string fmt;
            if (containment_res == ContainmentResult::Identical) {
              fmt = "Identical counted objects are not supported yet, error for %s and %s." ;
            }
            else {
              fmt = "Error while of counted object is not supported yet, error for %s and %s.";
            }
            mcell_warn(
                fmt.c_str(),
                counted_objects[obj1].get_geometry_object(world).name.c_str(),
                counted_objects[obj2].get_geometry_object(world).name.c_str()
            );
            res = false;
          }
          break;

        default:
          assert(false);
      }
    }
  }

  return res;
}


static const GeometryObject* get_direct_parent(
    const World* world, const GeomObjectInfo& obj_info, const ContainmentMap& contained_in_mapping) {

  auto it = contained_in_mapping.find(obj_info);
  if (it == contained_in_mapping.end()) {
    return nullptr;
  }
  const uint_set<GeomObjectInfo>& obj_contained_in = it->second;

  // need to find an object that is a direct parent
  //   p in contained_in_mapping(obj)
  // such that:
  //  contained_in_mapping(p) == contained_in_mapping(obj) - p
  //
  // i.e.: p is the direct intermediate between our object and all other objects obj is contained in
  //
  for (const GeomObjectInfo& parent: obj_contained_in) {

    // copy 'parents' and remove 'p'
    uint_set<GeomObjectInfo> obj_contained_in_less_p = obj_contained_in;
    obj_contained_in_less_p.erase_existing(parent);

    auto it = contained_in_mapping.find(parent);
    if (it == contained_in_mapping.end()) {
      // this object has no parent so it must be the direct parent
      return &parent.get_geometry_object(world);
    }
    const uint_set<GeomObjectInfo>& p_contained_in = it->second;

    if (obj_contained_in_less_p == p_contained_in) {
      // this is the closest parent
      return &parent.get_geometry_object(world);
    }
  }

  // nothing found - is outside of all objects
  assert(obj_contained_in.empty());
  return nullptr;
}


static counted_volume_index_t get_counted_volume_index_using_parents(
    World* world,
    const geometry_object_id_t obj_id,
    const ContainmentMap& contained_in_mapping,
    const IntersectingSet& intersecting_objects
) {
  if (intersecting_objects.count(obj_id) != 0) {
    return COUNTED_VOLUME_INDEX_INTERSECTS;
  }

  // search in this map uses only partition id
  GeomObjectInfo obj_info(PARTITION_ID_INVALID, obj_id);

  const uint_set<GeomObjectInfo>* obj_contained_in;

  auto it = contained_in_mapping.find(obj_info);
  if (it == contained_in_mapping.end()) {
    // not found - has no parent (intersection already checked)
    obj_contained_in = nullptr;
  }
  else {
    obj_contained_in = &it->second;
  }

  // we need to find or add a new counted volume to the partition
  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  // first prepare the counted volume
  CountedVolume cv;
  cv.contained_in_objects.insert(obj_id);
  if (obj_contained_in != nullptr) {
    for (const GeomObjectInfo& info: *obj_contained_in) {

      // if any of the parents intersects, let's assume that for counting this object intersects
      // as well, maybe some optimization can be done in the future
      if (intersecting_objects.count(info.geometry_object_id) != 0) {
        return COUNTED_VOLUME_INDEX_INTERSECTS;
      }

      geometry_object_index_t index = p.get_geometry_object(info.geometry_object_id).index;
      cv.contained_in_objects.insert_unique(index);
    }
  }

  return p.find_or_add_counted_volume(cv);
}


static void define_counted_volumes(
    World* world,
    const GeomObjectInfoVector& counted_objects,
    const ContainmentMap& contained_in_mapping,
    const IntersectingSet& intersecting_objects
) {
  for (const GeomObjectInfo& obj_info: counted_objects) {

    GeometryObject& current_obj = obj_info.get_geometry_object_noconst(world);

    // set inside index for our object
    if (intersecting_objects.count(obj_info.geometry_object_id) != 0) {
      // intersects with other objects - not sure in which counted volume
      // a molecule ends up when entering this object
      current_obj.counted_volume_index_inside = COUNTED_VOLUME_INDEX_INTERSECTS;
    }
    else {
      // does not intersect - get counted volume from its parents
      current_obj.counted_volume_index_inside =
          get_counted_volume_index_using_parents(
              world, obj_info.geometry_object_id, contained_in_mapping, intersecting_objects);
    }

    const GeometryObject* direct_parent_obj = get_direct_parent(world, obj_info, contained_in_mapping);

    // set outside index for our object
    counted_volume_index_t outside_index;
    // parent was found
    if (direct_parent_obj != nullptr) {
      outside_index =
          get_counted_volume_index_using_parents(
              world, direct_parent_obj->id, contained_in_mapping, intersecting_objects);
    }
    else {
      // parent is unclear
      if (intersecting_objects.count(obj_info.geometry_object_id) == 0) {
        // this is a top level object that has no parent but still needs a counted volume to be created
        outside_index = COUNTED_VOLUME_INDEX_OUTSIDE_ALL;
      }
      else {
        // it is not clear what the outside is
        outside_index = COUNTED_VOLUME_INDEX_INTERSECTS;
      }
    }

    current_obj.counted_volume_index_outside = outside_index;
  }
}


// return true if counted volumes were correctly set up
// the only entry point to this utility file (so far)
bool initialize_counted_volumes(World* world) {

  GeomObjectInfoVector counted_objects;

  // get counted objects
  for (const Partition& p: world->get_partitions()) {
    for (const GeometryObject& obj: p.get_geometry_objects()) {
      if (obj.is_counted_volume) {
        counted_objects.push_back( GeomObjectInfo(p.id, obj.id) );
      }
    }
  }

  // prepare VTK polydata objects
  convert_objects_to_clean_polydata(world, counted_objects);

  // holds mapping that tells in which geometry objects is a given object contained
  ContainmentMap contained_in_mapping;

  // set of object ids that intersect
  IntersectingSet intersecting_objects;

  // compute 'contained-in' mapping
  bool ok = compute_containement_mapping(world, counted_objects, contained_in_mapping, intersecting_objects);
  if (!ok) {
    return false;
  }

  // define counted volumes, sets inside and outside volume id for geometry objects
  uint_set<geometry_object_id_t> already_processed_objects;
  define_counted_volumes(world, counted_objects, contained_in_mapping, intersecting_objects);

  return true;
}


bool is_point_inside_counted_volume(GeometryObject& obj, const Vec3& point) {
  assert(obj.is_counted_volume);
  assert(obj.counted_volume_polydata.Get() != nullptr);

  double point_coords[3] = {point.x, point.y, point.z};

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(point_coords);

  vtkSmartPointer<vtkPolyData> points_polydata = vtkSmartPointer<vtkPolyData>::New();
  points_polydata->SetPoints(points);

  auto select_enclosed_points = vtkSmartPointer<vtkSelectEnclosedPoints>::New();

  select_enclosed_points->SetSurfaceData(obj.counted_volume_polydata);
  select_enclosed_points->SetInputData(points_polydata); // and we are trying whether poly1 fits into that
  select_enclosed_points->Update();

  return select_enclosed_points->IsInside(0) ;
}



} // namespace CountedVolumesUtil
} // namespace MCell
