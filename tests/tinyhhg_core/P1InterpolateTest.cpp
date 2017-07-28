#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include <tinyhhg_core/p1functionspace/p1face.hpp>
#include <tinyhhg_core/p1functionspace/p1edge.hpp>


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


  const size_t minLevel = 2;
  const size_t maxLevel = 5;

  size_t v_perFace = levelinfo::num_microvertices_per_face(maxLevel);
  size_t v_perEdge = levelinfo::num_microvertices_per_edge(maxLevel);
  size_t nbr_v_perEdge = v_perEdge - 1;
  size_t v_perVertex = levelinfo::num_microvertices_per_vertex(maxLevel);

  P1Function x("x", storage, minLevel, maxLevel);

  for (auto face : storage->getFaces())
  {
    for (size_t i = 0; i < v_perFace; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(face.second->getData(x.getFaceDataID())->data[maxLevel][i], 0.0);
    }
  }
  for (auto edge : storage->getEdges())
  {
    for (size_t i = 0; i < v_perEdge + edge.second->getNumHigherDimNeighbors() * nbr_v_perEdge; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(edge.second->getData(x.getEdgeDataID())->data[maxLevel][i], 0.0);
    }
  }
  for (auto vertex : storage->getVertices())
  {
    for (size_t i = 0; i < v_perVertex + vertex.second->getNumHigherDimNeighbors(); ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(vertex.second->getData(x.getVertexDataID())->data[maxLevel][i], 0.0);
    }
  }



  std::function<real_t(const Point3D &)> exact = [](const Point3D & xx) { return 2*xx[0] + xx[1]; };
  //std::function<real_t(const Point3D &)> exact = [](const Point3D & x) { return 13; };

  real_t value,xStepSize, yStepSize;
  for(auto faceIter : storage->getFaces())
  {
    Face *face = faceIter.second;
    value = face->coords[0].x[0] * 2 + face->coords[0].x[0];
//    Edge *faceEdge0 = storage->getEdge(face->getEdgeID0().getID());
//    Edge *faceEdge1 = storage->getEdge(face->getEdgeID1().getID());
    xStepSize = walberla::real_c(face->coords[1].x[0] - face->coords[0].x[0]) / walberla::real_c((v_perEdge-1));
    yStepSize = walberla::real_c(face->coords[2].x[1] - face->coords[0].x[1]) / walberla::real_c((v_perEdge-1));

    P1Face::interpolate(maxLevel, *face, x.getFaceDataID(), exact);
    for (uint_t i = 0; i < v_perEdge; ++i)
    {
      for (uint_t j = 0; j < v_perEdge - i; ++j)
      {
        uint_t idx = P1Face::index<maxLevel>(j, i, P1Face::C);
        if (P1Face::is_boundary(idx, v_perEdge))
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->data[maxLevel][idx], 0.0,
                                     "i: " << i << " j: " << j << " idx: " << idx << " value was " << value);
        } else
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->data[maxLevel][idx], value,
                                     "i: " << i << " j: " << j << " idx: " << idx << " value was " << value);
        }
        value += 2 * xStepSize;
      }
      value = yStepSize * (i + 1);
    }
  }

  value = 0;
  for(auto edgeIter : storage->getEdges()){
    Edge* edge = edgeIter.second;
    hhg::P1Edge::interpolate(*edge,x.getEdgeDataID(),exact,maxLevel);
    value = 2 * edge->coords_[0].x[0] + edge->coords_[0].x[1];
    xStepSize = edge->direction_.x[0] / walberla::real_c((v_perEdge-1));
    yStepSize = edge->direction_.x[1] / walberla::real_c((v_perEdge-1));
    value += 2*xStepSize;
    value += yStepSize;
    for(uint_t i = 1; i < v_perEdge-1; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(edge->getData(x.getEdgeDataID())->data[maxLevel][i], value,
                                 "i: " << i << " edge: "<< *edge);
      value += 2*xStepSize;
      value += yStepSize;
    }
  }

  return EXIT_SUCCESS;
}
