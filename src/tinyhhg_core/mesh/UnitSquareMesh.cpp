#include "MeshInfo.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

namespace hhg {

MeshInfo MeshInfo::unitSquareMesh(uint_t level) {
  MeshInfo meshInfo;

  uint_t N = (uint_t) std::pow(2 ,level);
  real_t h = walberla::real_c(1.0)/walberla::real_c(N);

  for (uint_t row = 0; row < N+1; ++row) {
    for (uint_t col = 0; col < N+1; ++col) {
      uint_t id = (N+1)*row + col;
      meshInfo.vertices_[id] =  MeshInfo::Vertex(id, h * col * Point3D{{{1.0, 0.0, 0.0}}} + h * row * Point3D{{{0.0, 1.0, 0.0}}}, Inner);
    }
  }

  std::vector<IDType> triangleVertices(3);

  for (uint_t row = 0; row < N; ++row) {
    for (uint_t col = 0; col < N; ++col) {

      // Up triangle
      triangleVertices[0] = (N+1) * row + col;
      triangleVertices[1] = (N+1) * row + (col+1);
      triangleVertices[2] = (N+1) * (row+1) + col;

      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({{triangleVertices[0], triangleVertices[1]}}), (row == 0) ? DirichletBoundary : Inner));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({{triangleVertices[1], triangleVertices[2]}}), Inner));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({{triangleVertices[2], triangleVertices[0]}}), (col == 0) ? DirichletBoundary : Inner));

      meshInfo.addFace(MeshInfo::Face(triangleVertices, Inner));

      // Down  triangle
      triangleVertices[0] = (N+1) * (row+1) + (col+1);
      triangleVertices[1] = (N+1) * (row+1) + col;
      triangleVertices[2] = (N+1) * row + (col+1);

      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({{triangleVertices[0], triangleVertices[1]}}), (row == N-1) ? DirichletBoundary : Inner));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({{triangleVertices[1], triangleVertices[2]}}), Inner));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({{triangleVertices[2], triangleVertices[0]}}), (col == N-1) ? DirichletBoundary : Inner));

      meshInfo.addFace(MeshInfo::Face(triangleVertices, Inner));
    }
  }

  return meshInfo;
}

}