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

#include "data_model_geometry.h"

#include <stdexcept>

#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkTransform.h>
#include <vtkPointData.h>
#include <vtkFeatureEdges.h>
#include <vtkMassProperties.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>

#include "include/datamodel_defines.h"
#include "generator_utils.h"

using namespace Json;
using Json::Value;

namespace MCell {

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


static vtkSmartPointer<vtkPolyData> convert_dm_object_to_polydata(
    Json::Value& model_object) {

  string name = model_object[KEY_NAME].asString();

  // we need to convert each geometry object into VTK's polydata representation
  // example: https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/TriangleArea
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

  Value& vertex_list = get_node(name, model_object, KEY_VERTEX_LIST);
  Value& element_connections = get_node(name, model_object, KEY_ELEMENT_CONNECTIONS);

  // assuming the the data model is correct
  for (Value::ArrayIndex i = 0; i < vertex_list.size(); i++) {
    Value& vertex = vertex_list[i];
    points->InsertNextPoint(vertex[0].asDouble(), vertex[1].asDouble(), vertex[2].asDouble());
  }

  // store triangles
  for (Value::ArrayIndex i = 0; i < element_connections.size(); i++) {
    Value& element = element_connections[i];
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      triangle->GetPointIds()->SetId(i, element[i].asInt());
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


// practically the same as in vtk_utils
// auxiliary function to compute volume, not related to counted volumes but uses VTK
// returns FLT_INVALID if the object is not watertight
static void get_polydata_volume_and_area(
    vtkSmartPointer<vtkPolyData> polydata, double& volume, double& area) {

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

  volume = massProperties->GetVolume();
  area = massProperties->GetSurfaceArea();
}


void compute_volume_and_area(Json::Value& model_object, double& volume, double& area) {

  string name = model_object[KEY_NAME].asString();

  vtkSmartPointer<vtkPolyData> polydata = convert_dm_object_to_polydata(model_object);

  if (!is_watertight(polydata)) {
    throw ConversionError("Geometry object " + name + " is not watertight, could not compute volume.");
  }

  get_polydata_volume_and_area(polydata, volume, area);
}

} // namespace MCell
