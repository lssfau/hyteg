#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <core/Environment.h>
#include <core/config/Config.h>

#include <tinyhhg_core/geometry/Geometry.hpp>
#include <tinyhhg_core/geometry/IdentityMap.hpp>
#include <tinyhhg_core/geometry/CircularMap.hpp>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{

  /// create enviroment
  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  const size_t level = 5;

  /// read mesh file and create storage
  MeshInfo meshInfo = MeshInfo::fromGmshFile( "../data/meshes/unitsquare_with_circular_hole.msh" );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  hhg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> test = [](const hhg::Point3D&,const std::vector<real_t>&) { return 1.0; };

  auto x = std::make_shared<hhg::P1Function< real_t > >("x", storage, level, level);
  auto y = std::make_shared<hhg::P1Function< real_t > >("y", storage, level, level);

  auto k11 = std::make_shared<hhg::P1Function< real_t > >("k11", storage, level, level);
  auto k12 = std::make_shared<hhg::P1Function< real_t > >("k12", storage, level, level);
  auto k22 = std::make_shared<hhg::P1Function< real_t > >("k22", storage, level, level);

  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Vertex>> emptyVertexIds;
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Edge>> emptyEdgeIds;
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Face>> emptyFaceIds;

  hhg::Geometry geometry(storage);

  for (auto& it : storage->getFaces()) {
    Face& face = *it.second;

    auto geometryMap = face.getData(geometry.getFaceMapID());

    // Create map
    if (face.hasBoundaryEdge()) {
      WALBERLA_LOG_INFO("Detected boundary face");
      geometryMap->map = std::shared_ptr<FaceMap2D>(new CircularMap(face, storage, {{0.5, 0.5}}, 0.25));
    } else {
      geometryMap->map = std::shared_ptr<FaceMap2D>(new IdentityMap(face));
    }

    // Apply map
    std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> tmp_x = [&](const hhg::Point3D& x_,const std::vector<real_t>&) {
      hhg::Point2D tmp, out;
      tmp[0] = x_[0];
      tmp[1] = x_[1];
      geometryMap->map->evalF(tmp, out);
      return out[0];
    };

    std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> tmp_y = [&](const hhg::Point3D& x_,const std::vector<real_t>&) {
      hhg::Point2D tmp, out;
      tmp[0] = x_[0];
      tmp[1] = x_[1];
      geometryMap->map->evalF(tmp, out);
      return out[1];
    };

    std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> tmp_k11 = [&](const hhg::Point3D& x_,const std::vector<real_t>&) {
      hhg::Point2D tmp;
      hhg::Matrix2r out;
      tmp[0] = x_[0];
      tmp[1] = x_[1];
      geometryMap->map->evalDF(tmp, out);
      return out(0,0);
    };

    std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> tmp_k12 = [&](const hhg::Point3D& x_,const std::vector<real_t>&) {
      hhg::Point2D tmp;
      hhg::Matrix2r out;
      tmp[0] = x_[0];
      tmp[1] = x_[1];
      geometryMap->map->evalDF(tmp, out);
      return out(0,1);
    };

    std::function<real_t(const hhg::Point3D&,const std::vector<real_t>&)> tmp_k22 = [&](const hhg::Point3D& x_,const std::vector<real_t>&) {
      hhg::Point2D tmp;
      hhg::Matrix2r out;
      tmp[0] = x_[0];
      tmp[1] = x_[1];
      geometryMap->map->evalDF(tmp, out);
      return out(1,1);
    };

    vertexdof::macroface::interpolate< real_t >(level, face, x->getFaceDataID(), emptyFaceIds, tmp_x);
    vertexdof::macroface::interpolate< real_t >(level, face, y->getFaceDataID(), emptyFaceIds, tmp_y);
    vertexdof::macroface::interpolate< real_t >(level, face, k11->getFaceDataID(), emptyFaceIds, tmp_k11);
    vertexdof::macroface::interpolate< real_t >(level, face, k12->getFaceDataID(), emptyFaceIds, tmp_k12);
    vertexdof::macroface::interpolate< real_t >(level, face, k22->getFaceDataID(), emptyFaceIds, tmp_k22);

    for (auto& edgeId : face.neighborEdges()) {
      Edge& edge = *storage->getEdge(edgeId);
      vertexdof::macroedge::interpolate< real_t >(level, edge, x->getEdgeDataID(), emptyEdgeIds, tmp_x);
      vertexdof::macroedge::interpolate< real_t >(level, edge, y->getEdgeDataID(), emptyEdgeIds, tmp_y);
      vertexdof::macroedge::interpolate< real_t >(level, edge, k11->getEdgeDataID(), emptyEdgeIds, tmp_k11);
      vertexdof::macroedge::interpolate< real_t >(level, edge, k12->getEdgeDataID(), emptyEdgeIds, tmp_k12);
      vertexdof::macroedge::interpolate< real_t >(level, edge, k22->getEdgeDataID(), emptyEdgeIds, tmp_k22);
    }

    for (auto& vertexId : face.neighborVertices()) {
      Vertex& vertex = *storage->getVertex(vertexId);
      vertexdof::macrovertex::interpolate< real_t >(vertex, x->getVertexDataID(), emptyVertexIds, tmp_x, level);
      vertexdof::macrovertex::interpolate< real_t >(vertex, y->getVertexDataID(), emptyVertexIds, tmp_y, level);
      vertexdof::macrovertex::interpolate< real_t >(vertex, k11->getVertexDataID(), emptyVertexIds, tmp_k11, level);
      vertexdof::macrovertex::interpolate< real_t >(vertex, k12->getVertexDataID(), emptyVertexIds, tmp_k12, level);
      vertexdof::macrovertex::interpolate< real_t >(vertex, k22->getVertexDataID(), emptyVertexIds, tmp_k22, level);
    }
  }

  VTKOutput vtkOutput("../output", "GeometryBlending");
  vtkOutput.add(x.get());
  vtkOutput.add(y.get());
  vtkOutput.add(k11.get());
  vtkOutput.add(k12.get());
  vtkOutput.add(k22.get());
  vtkOutput.write(level);

  return 0;
}
