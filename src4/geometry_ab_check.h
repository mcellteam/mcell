/******************************************************************************
 *
 * Copyright (C) 2026 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

// Temporary, for the VTK removal: built with MCELL_GEOMETRY_AB_CHECK, every
// place that uses VTK also runs geometry_analysis on the same mesh and records
// the comparison. A difference is printed to stderr as one line starting with
// "MCELL_GEOMETRY_AB:"; when the environment variable MCELL_GEOMETRY_AB_LOG_DIR
// is set, every comparison is also appended to <dir>/ab_<pid>.log, so that a
// test run can show that the checks ran. Removed together with VTK.

#ifndef SRC4_GEOMETRY_AB_CHECK_H_
#define SRC4_GEOMETRY_AB_CHECK_H_

#ifdef MCELL_GEOMETRY_AB_CHECK

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

#include <vtkCellArray.h>
#include <vtkPolyData.h>

#include "geometry_analysis.h"

namespace MCell {
namespace GeometryAbCheck {

inline long process_id() {
#ifdef _WIN32
  return _getpid();
#else
  return getpid();
#endif
}

// an empty difference means that VTK and geometry_analysis agree
inline void record(const std::string& site, const std::string& what, const std::string& difference) {
  if (!difference.empty()) {
    fprintf(stderr, "MCELL_GEOMETRY_AB: %s: %s: %s\n", site.c_str(), what.c_str(), difference.c_str());
    fflush(stderr);
  }
  const char* dir = std::getenv("MCELL_GEOMETRY_AB_LOG_DIR");
  if (dir != nullptr && dir[0] != '\0') {
    const std::string path = std::string(dir) + "/ab_" + std::to_string(process_id()) + ".log";
    FILE* f = fopen(path.c_str(), "a");
    if (f != nullptr) {
      fprintf(f, "%s\t%s\t%s\n", site.c_str(), what.c_str(), difference.empty() ? "same" : ("DIFF " + difference).c_str());
      fclose(f);
    }
  }
}

inline std::string compare(bool vtk, bool ga) {
  return vtk == ga ? "" : std::string("vtk ") + (vtk ? "true" : "false") + ", geometry_analysis " + (ga ? "true" : "false");
}

// compared bit for bit
inline std::string compare(double vtk, double ga) {
  if (std::memcmp(&vtk, &ga, sizeof(double)) == 0) {
    return "";
  }
  char buf[128];
  snprintf(buf, sizeof(buf), "vtk %.17g, geometry_analysis %.17g", vtk, ga);
  return buf;
}

// VTK's cleaned polydata against the module's cleaned mesh: the same points
// in the same order and the same triangles
inline std::string compare(vtkPolyData* pd, const GeometryAnalysis::Mesh& mesh) {
  if (pd->GetNumberOfLines() != 0 || pd->GetNumberOfVerts() != 0 || pd->GetNumberOfStrips() != 0) {
    return "VTK turned degenerate triangles into lines or vertices";
  }
  if (pd->GetNumberOfPoints() != static_cast<vtkIdType>(mesh.points.size())) {
    return "vtk " + std::to_string(pd->GetNumberOfPoints()) + " points, geometry_analysis " +
        std::to_string(mesh.points.size());
  }
  for (vtkIdType i = 0; i < pd->GetNumberOfPoints(); i++) {
    double p[3];
    pd->GetPoint(i, p);
    if (p[0] != mesh.points[i][0] || p[1] != mesh.points[i][1] || p[2] != mesh.points[i][2]) {
      return "point " + std::to_string(i) + " differs";
    }
  }
  vtkCellArray* polys = pd->GetPolys();
  if (polys->GetNumberOfCells() != static_cast<vtkIdType>(mesh.triangles.size())) {
    return "vtk " + std::to_string(polys->GetNumberOfCells()) + " triangles, geometry_analysis " +
        std::to_string(mesh.triangles.size());
  }
  vtkIdType npts;
  const vtkIdType* pts;
  size_t t = 0;
  for (polys->InitTraversal(); polys->GetNextCell(npts, pts); t++) {
    if (npts != 3 || pts[0] != mesh.triangles[t][0] || pts[1] != mesh.triangles[t][1] ||
        pts[2] != mesh.triangles[t][2]) {
      return "triangle " + std::to_string(t) + " differs";
    }
  }
  return "";
}

} // namespace GeometryAbCheck
} // namespace MCell

#endif // MCELL_GEOMETRY_AB_CHECK

#endif // SRC4_GEOMETRY_AB_CHECK_H_
