#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../../data/meshes/tri_1el_neumann.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(meshFileName);
  hhg::SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t level = 2;
  size_t maxiter = 10000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::P2Function<real_t> u("u", storage, level, level);
  std::function<real_t(const hhg::Point3D&)> testExpression = [](const hhg::Point3D& x) { return x[1] + 1.0; };
  u.interpolate(testExpression, level, hhg::All);

  // Sync interpolated function values
  u.getCommunicator( level )->template communicate< Vertex, Edge >();
  u.getCommunicator( level )->template communicate< Edge, Face >();
  u.getCommunicator( level )->template communicate< Face, Edge >();
  u.getCommunicator( level )->template communicate< Edge, Vertex >();

  // We assume that the bottom left vertex has following connectivity
  //
  // EdgeDoF(1)
  //    |
  //    |   EdgeDoF(2)
  //    |   /
  // Vertex(0) --- EdgeDoF(0)

  // Get bottom left vertex
  auto vertex = storage->getVertex(PrimitiveID(0));

  // Get vertex values
  auto vertexEdgeData = vertex->getData(u.getEdgeDoFFunction()->getVertexDataID())->getPointer( level );

  WALBERLA_CHECK_EQUAL( vertexEdgeData[0], 1.0 );
  WALBERLA_CHECK_EQUAL( vertexEdgeData[1], 1.125 );
  WALBERLA_CHECK_EQUAL( vertexEdgeData[2], 1.125 );

  return EXIT_SUCCESS;
}
