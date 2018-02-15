#include <core/timing/Timer.h>
#include <tinyhhg_core/tinyhhg.hpp>
#include <core/Environment.h>
#include <core/config/Config.h>

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

  std::function<real_t(const hhg::Point3D&)> test = [](const hhg::Point3D&) { return 1.0; };

  auto x = std::make_shared<hhg::P1Function< real_t > >("x", storage, level, level);
  auto y = std::make_shared<hhg::P1Function< real_t > >("y", storage, level, level);

  std::function<real_t(const hhg::Point3D&)> tmp_x = [&](const hhg::Point3D& x_) {
    return x_[0];
  };

  std::function<real_t(const hhg::Point3D&)> tmp_y = [&](const hhg::Point3D& x_) {
    return x_[1];
  };

  Point3D circleCenter{{0.5, 0.5, 0.0}};
  real_t circleRadius = 0.25;

  for (auto& it : storage->getFaces()) {
    Face &face = *it.second;

    if (face.hasBoundaryEdge()) {

      for (auto edgeId : face.edgesOnBoundary) {
        Edge& edge = *storage->getEdge(edgeId);

        if ((edge.getCoordinates()[0] - circleCenter).norm() < 0.9) {
          edge.setBlendingMap(std::shared_ptr<FaceMap>(new CircularMap(face, storage, circleCenter, circleRadius)));
        }
      }

      face.setBlendingMap(std::shared_ptr<FaceMap>(new CircularMap(face, storage, circleCenter, circleRadius)));
    }
  }

  x->interpolate(tmp_x, level, hhg::All);
  y->interpolate(tmp_y, level, hhg::All);

  VTKOutput vtkOutput("../output", "GeometryBlending");
  vtkOutput.add(x.get());
  vtkOutput.add(y.get());
  vtkOutput.write(level);

  return 0;
}
