#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include <tinyhhg_core/p1functionspace/p1face.hpp>
#include <tinyhhg_core/p1functionspace/p1edge.hpp>

#include <memory>

using walberla::real_t;
using namespace hhg;

int main(int argc, char **argv)
{
  walberla::debug::enterTestMode();
  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);


  size_t minLevel = 2;
  size_t maxLevel = 5;

  size_t v_perFace = levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = levelinfo::num_microvertices_per_edge(maxLevel);
  size_t nbr_v_perEdge = v_perEdge - 1;
  size_t v_perVertex = levelinfo::num_microvertices_per_vertex(maxLevel);

  P1Function x("x", storage, minLevel, maxLevel);

  for (auto face : storage.get()->getFaces())
  {
    for (size_t i = 0; i < v_perFace; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(face.second->getData(x.getFaceDataID())->data[maxLevel][i], 0.0);
    }
  }
  for (auto edge : storage.get()->getEdges())
  {
    for (size_t i = 0; i < v_perEdge + edge.second->getNumHigherDimNeighbors() * nbr_v_perEdge; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(edge.second->getData(x.getEdgeDataID())->data[maxLevel][i], 0.0);
    }
  }
  for (auto vertex : storage.get()->getVertices())
  {
    for (size_t i = 0; i < v_perVertex + vertex.second->getNumHigherDimNeighbors(); ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(vertex.second->getData(x.getVertexDataID())->data[maxLevel][i], 0.0);
    }
  }

  //since our triangle is rectangular and at (0,0),(1,0) and (0,1) we can calculate easily
  real_t xStepSize = 1 / (v_perEdge-1);
  real_t yStepSize = 1 / (v_perEdge-1);

  std::function<real_t(const Point3D &)> exact = [](const Point3D & x) { return x[0] + x[1]; };

  Face *faceZero = (*storage.get()->beginFaces()).second;
  P1Face::interpolate(maxLevel,*faceZero, x.getFaceDataID(), exact);


  for (size_t i = 0; i < v_perFace; ++i)
  {
    if (P1Face::is_boundary(i, v_perEdge))
    {
      WALBERLA_CHECK_FLOAT_EQUAL(faceZero->getData(x.getFaceDataID())->data[maxLevel][i], 0.0);
    } else
    {
      WALBERLA_CHECK_FLOAT_EQUAL(faceZero->getData(x.getFaceDataID())->data[maxLevel][i], 13.0);
    }
  }

  Edge *edgeZero = (*storage.get()->beginEdges()).second;
  P1Edge::interpolate(*edgeZero, x.getEdgeDataID(), exact, maxLevel);

  for (size_t i = 1; i < (v_perEdge - 1); ++i)
  {
    WALBERLA_CHECK_FLOAT_EQUAL(edgeZero->getData(x.getEdgeDataID())->data[maxLevel][i], 13.0)
  }
  WALBERLA_CHECK_FLOAT_EQUAL(edgeZero->getData(x.getEdgeDataID())->data[maxLevel][0], 0.0)
  WALBERLA_CHECK_FLOAT_EQUAL(edgeZero->getData(x.getEdgeDataID())->data[maxLevel][v_perEdge - 1], 0.0)

  for (size_t i = v_perEdge; i < v_perEdge + edgeZero->getNumHigherDimNeighbors() * nbr_v_perEdge; ++i)
  {
    WALBERLA_CHECK_FLOAT_EQUAL(edgeZero->getData(x.getEdgeDataID())->data[maxLevel][i], 0.0);
  }

  return EXIT_SUCCESS;
}
