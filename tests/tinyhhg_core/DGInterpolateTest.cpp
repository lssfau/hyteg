#include "tinyhhg_core/dgfunctionspace/DGFunction.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

#include "core/Environment.h"

using walberla::real_t;
using namespace hhg;

int main(int argc, char **argv) {
  walberla::debug::enterTestMode();
  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);


  const size_t minLevel = 2;
  const size_t maxLevel = 5;

  size_t v_perFace = levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = levelinfo::num_microvertices_per_edge(maxLevel);
  size_t nbr_v_perEdge = v_perEdge - 1;
  size_t v_perVertex = levelinfo::num_microvertices_per_vertex(maxLevel);

  std::function<real_t(const Point3D &)> exact = [](const Point3D & xx) { return 2*xx[0] + xx[1]; };

  std::shared_ptr<Edge> edge0 = (*storage->getEdges().begin()).second;
  DGFunction < real_t > x("x", storage, minLevel, maxLevel);

  DGEdge::interpolate< real_t >(maxLevel,*edge0,x.getEdgeDataID(),exact,storage);

  x.interpolate(exact,maxLevel);

}