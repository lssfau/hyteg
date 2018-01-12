
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"
#include "core/debug/all.h"

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
  size_t mEperEdge = levelinfo::num_microedges_per_edge(maxLevel);
  size_t nbr_v_perEdge = mEperEdge - 1;
  size_t v_perVertex = levelinfo::num_microvertices_per_vertex(maxLevel);

  EdgeDoFFunction< real_t > x("x", storage, minLevel, maxLevel);
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Edge>> emptyEdgeIds;
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Face>> emptyFaceIds;

//  for (auto face : storage->getFaces())
//  {
//    for (size_t i = 0; i < v_perFace; ++i)
//    {
//      WALBERLA_CHECK_FLOAT_EQUAL(face.second->getData(x.getFaceDataID())->getPointer(maxLevel)[i], 0.0);
//    }
//  }
//  for (auto edge : storage->getEdges())
//  {
//    for (size_t i = 0; i < mEperEdge + edge.second->getNumHigherDimNeighbors() * nbr_v_perEdge; ++i)
//    {
//      WALBERLA_CHECK_FLOAT_EQUAL(edge.second->getData(x.getEdgeDataID())->getPointer(maxLevel)[i], 0.0);
//    }
//  }
//  for (auto vertex : storage->getVertices())
//  {
//    for (size_t i = 0; i < v_perVertex + vertex.second->getNumHigherDimNeighbors(); ++i)
//    {
//      WALBERLA_CHECK_FLOAT_EQUAL(vertex.second->getData(x.getVertexDataID())->getPointer(maxLevel)[i], 0.0);
//    }
//  }



  std::function<real_t(const Point3D &,const std::vector<real_t>&)> exact = [](const Point3D & xx, const std::vector<real_t>&) { return 2*xx[0] + xx[1]; };
  //std::function<real_t(const Point3D &)> exact = [](const Point3D & x) { return 13; };

  real_t value,xStepSize, yStepSize;
  for(auto faceIter : storage->getFaces())
  {
    auto face = faceIter.second;

//    Edge *faceEdge0 = storage->getEdge(face->getEdgeID0().getID());
//    Edge *faceEdge1 = storage->getEdge(face->getEdgeID1().getID());
    xStepSize = walberla::real_c(face->coords[1].x[0] - face->coords[0].x[0]) / walberla::real_c((mEperEdge));
    yStepSize = walberla::real_c(face->coords[2].x[1] - face->coords[0].x[1]) / walberla::real_c((mEperEdge));
    value = (xStepSize / 2 + face->coords[0].x[0]) * 2 + face->coords[0].x[0];

    edgedof::macroface::interpolate< real_t >(maxLevel, *face, x.getFaceDataID(), emptyFaceIds, exact);
    for (uint_t i = 0; i < mEperEdge; ++i)
    {
      for (uint_t j = 0; j < mEperEdge - i; ++j)
      {
        uint_t idx_ho = indexing::edgedof::macroface::indexFromHorizontalEdge<maxLevel>(j, i, stencilDirection::EDGE_HO_C);
        uint_t idx_di = indexing::edgedof::macroface::indexFromDiagonalEdge<maxLevel>(j, i, stencilDirection::EDGE_DI_C);
        uint_t idx_ve = indexing::edgedof::macroface::indexFromVerticalEdge<maxLevel>(j, i, stencilDirection::EDGE_VE_C);
        if (i == 0)
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx_ho], 0.0,
                                     "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value);
        } else
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx_ho], value,
                                     "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value);
        }
        if (j == 0)
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx_ve], 0.0,
                                     "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value);
        } else
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx_ve], value + yStepSize / 2 - xStepSize,
                                     "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value + yStepSize / 2 - xStepSize);
        }
        if (i + j == levelinfo::num_microvertices_per_edge(maxLevel) - 2)
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx_di], 0.0,
                                     "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value);
        } else
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx_di], value + yStepSize / 2 ,
                                     "i: " << i << " j: " << j << " idx: " << idx_ho << " value was " << value + yStepSize / 2);
        }
        value += 2 * xStepSize;
      }
      value = yStepSize * (real_t) (i + 1) + xStepSize / 2 * 2;
    }
  }

  value = 0;
  for(auto edgeIter : storage->getEdges()){
    auto edge = edgeIter.second;
    hhg::edgedof::macroedge::interpolate< real_t >(maxLevel, *edge,x.getEdgeDataID(),emptyEdgeIds,exact);
    value = 2 * edge->getCoordinates()[0].x[0] + edge->getCoordinates()[0].x[1];
    xStepSize = edge->getDirection().x[0] / walberla::real_c((mEperEdge-1));
    yStepSize = edge->getDirection().x[1] / walberla::real_c((mEperEdge-1));
    value += 2*xStepSize;
    value += yStepSize;
    for(uint_t i = 1; i < mEperEdge-1; ++i)
    {
//      WALBERLA_CHECK_FLOAT_EQUAL(edge->getData(x.getEdgeDataID())->getPointer(maxLevel)[i], value,
//                                 "i: " << i << " edge: "<< *edge);
      value += 2*xStepSize;
      value += yStepSize;
    }
  }

  return EXIT_SUCCESS;
}
