/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/


#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkCollisionDetectionFilter.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkTransform.h>
#include <vtkPointData.h>
#include <vtkFeatureEdges.h>
#include <vtkMassProperties.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkOBJExporter.h>

#include "logging.h"

#include "world.h"
#include "partition.h"
#include "geometry.h"

#include "vtk_utils.h"

using namespace std;

namespace MCell {

namespace VtkUtils {

enum class ContainmentResult {
  Error,
  Disjoint,
  Identical,
  Intersect,
  Obj1InObj2,
  Obj2InObj1
};


GeometryObject& GeomObjectInfo::get_geometry_object_noconst(World* world) const {
  assert(world != nullptr);
  assert(for_counted_objects);
  Partition& p = world->get_partition(partition_id);
  return p.get_geometry_object(geometry_object_id);
}


const GeometryObject& GeomObjectInfo::get_geometry_object(const World* world) const {
  assert(world != nullptr);
  assert(for_counted_objects);
  const Partition& p = world->get_partition(partition_id);
  return p.get_geometry_object(geometry_object_id);
}


static vtkSmartPointer<vtkPolyData> convert_geometry_object_to_polydata(
    const World* world, const GeometryObject& obj, const bool convert_units = false) {

  const Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  // we need to convert each geometry object into VTK's polydata representation
  // example: https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/TriangleArea
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

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
    if (convert_units) {
      pt = pt * Vec3(world->config.length_unit);
    }
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

  // also copy this information to our object
  return clean->GetOutput();
}


static bool convert_objects_to_clean_polydata(World* world, GeomObjectInfoVector& counted_objects) {
  bool res = true;

  for (GeomObjectInfo& obj_info: counted_objects) {

    GeometryObject& obj = world->get_geometry_object(obj_info.geometry_object_id);
    vtkSmartPointer<vtkPolyData> polydata = convert_geometry_object_to_polydata(world, obj);


    int closed = vtkSelectEnclosedPoints::IsSurfaceClosed(polydata.Get());
    if (closed != 1) {
      mcell_warn("Counting object must be closed, error for %s.", obj.name.c_str());
      res = false;
      continue;
    }

    // copy a reference to object info as well
    obj_info.polydata = polydata;
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
bool compute_containement_mapping(
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
          intersecting_objects.insert(counted_objects[obj1]);
          intersecting_objects.insert(counted_objects[obj2]);
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
            if (world != nullptr) {
              mcell_warn(
                  fmt.c_str(),
                  counted_objects[obj1].get_geometry_object(world).name.c_str(),
                  counted_objects[obj2].get_geometry_object(world).name.c_str()
              );
            }
            else {
              mcell_warn(
                  fmt.c_str(),
                  counted_objects[obj1].name.c_str(),
                  counted_objects[obj2].name.c_str()
              );
            }
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


const GeomObjectInfo* get_direct_parent_info(
    const GeomObjectInfo& obj_info, const ContainmentMap& contained_in_mapping) {

  auto it = contained_in_mapping.find(obj_info);
  if (it == contained_in_mapping.end()) {
    return nullptr;
  }
  const std::set<GeomObjectInfo>& obj_contained_in = it->second;

  // need to find an object that is a direct parent
  //   p in contained_in_mapping(obj)
  // such that:
  //  contained_in_mapping(p) == contained_in_mapping(obj) - p
  //
  // i.e.: p is the direct intermediate between our object and all other objects obj is contained in
  //
  for (const GeomObjectInfo& parent: obj_contained_in) {

    // copy 'parents' and remove 'p'
    std::set<GeomObjectInfo> obj_contained_in_less_p = obj_contained_in;
    obj_contained_in_less_p.erase(parent);

    auto it = contained_in_mapping.find(parent);
    if (it == contained_in_mapping.end()) {
      // this object has no parent so it must be the direct parent
      return &parent;
    }
    const std::set<GeomObjectInfo>& p_contained_in = it->second;

    if (obj_contained_in_less_p == p_contained_in) {
      // this is the closest parent
      return &parent;
    }
  }

  // nothing found - is outside of all objects
  assert(obj_contained_in.empty());
  return nullptr;
}


static const GeometryObject* get_direct_parent(
    const World* world, const GeomObjectInfo& obj_info, const ContainmentMap& contained_in_mapping) {

  const GeomObjectInfo* info = get_direct_parent_info(obj_info, contained_in_mapping);
  if (info == nullptr) {
    return nullptr;
  }
  else {
    return &info->get_geometry_object(world);
  }
}


static counted_volume_index_t get_counted_volume_index_using_parents(
    World* world,
    const geometry_object_id_t obj_id,
    const ContainmentMap& contained_in_mapping,
    const IntersectingSet& intersecting_objects
) {
  // search in this map uses only object id
  GeomObjectInfo obj_info(PARTITION_ID_INVALID, obj_id);

  if (intersecting_objects.count(obj_info) != 0) {
    return COUNTED_VOLUME_INDEX_INTERSECTS;
  }

  const std::set<GeomObjectInfo>* obj_contained_in;

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
      if (intersecting_objects.count(info) != 0) {
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
    if (intersecting_objects.count(obj_info) != 0) {
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
      if (intersecting_objects.count(obj_info) == 0) {
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
bool initialize_counted_volumes(World* world, bool& has_intersecting_counted_objects) {

  GeomObjectInfoVector counted_objects;

  // get counted objects
  for (const Partition& p: world->get_partitions()) {
    for (const GeometryObject& obj: p.get_geometry_objects()) {
      if (obj.is_counted_volume_or_compartment()) {
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
  define_counted_volumes(world, counted_objects, contained_in_mapping, intersecting_objects);

  has_intersecting_counted_objects = !intersecting_objects.empty();

  return true;
}

#if 0
bool is_point_inside_counted_volume(GeometryObject& obj, const Vec3& point) {
  assert(obj.is_counted_volume_or_compartment());
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
#endif


static double is_watertight(vtkSmartPointer<vtkPolyData> polydata) {
  // check if the object is watertight,
  // based on https://lorensen.github.io/VTKExamples/site/Cxx/Meshes/BoundaryEdges/
  vtkSmartPointer<vtkFeatureEdges> featureEdges =
    vtkSmartPointer<vtkFeatureEdges>::New();
  featureEdges->SetInputData(polydata.Get());
  featureEdges->BoundaryEdgesOn();
  featureEdges->FeatureEdgesOff();
  featureEdges->ManifoldEdgesOff();
  featureEdges->NonManifoldEdgesOff();
  featureEdges->Update();

  // if there are no edges, the object is watertight
  return featureEdges->GetOutput()->GetNumberOfCells() < 1;
}

// auxiliary function to compute volume, not related to counted volumes but uses VTK
// returns FLT_INVALID if the object is not watertight
double get_geometry_object_volume(const World* world, const GeometryObject& obj) {
  vtkSmartPointer<vtkPolyData> polydata = convert_geometry_object_to_polydata(world, obj);

  if (!is_watertight(polydata)) {
    return FLT_INVALID;
  }

  // volume computation is based on
  // https://lorensen.github.io/VTKExamples/site/Cxx/Utilities/MassProperties/
  vtkSmartPointer<vtkTriangleFilter> triangleFilter =
      vtkSmartPointer<vtkTriangleFilter>::New();
  triangleFilter->SetInputData(polydata);

  // Make the triangle windong order consistent
  vtkSmartPointer<vtkPolyDataNormals> normals =
    vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetInputConnection(triangleFilter->GetOutputPort());
  normals->ConsistencyOn();
  normals->SplittingOff();

  vtkSmartPointer<vtkMassProperties> massProperties =
    vtkSmartPointer<vtkMassProperties>::New();
  massProperties->SetInputConnection(normals->GetOutputPort());
  massProperties->Update();

  return massProperties->GetVolume();
}


void export_geometry_objects_to_obj(
    const World* world, const GeometryObjectVector& objs, const std::string& file_prefix) {

  vtkNew<vtkRenderer> renderer;

  vector<string> names;
  for (const GeometryObject& o: objs) {
    // convert objects
    vtkSmartPointer<vtkPolyData> polydata =
        convert_geometry_object_to_polydata(world, o, true);

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(polydata);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);

    // material, the default is white with 25% transparency
    // note: same as DEFAULT_COLOR
    actor->GetProperty()->SetDiffuseColor(1, 1, 1);
    actor->GetProperty()->SetOpacity(0.25);

    renderer->AddActor(actor);

    // we must also set names to the exported objects,
    // there is no support for object names in VTK and
    // the OBJ exporter needed to be updated
    names.push_back(o.name);
  }

  vtkNew<vtkRenderWindow> window;
  window->AddRenderer(renderer);

  vtkNew<vtkOBJExporter> exporter;
  exporter->SetRenderWindow(window);
  exporter->SetFilePrefix(file_prefix.c_str());
  exporter->SetActorNames(names); // custom MCell method
  exporter->Write();
}


} // namespace VtkUtils
} // namespace MCell
