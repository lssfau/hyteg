/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include <hyteg/p1functionspace/VertexDoFMacroEdge.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>
#include "core/Environment.h"

using walberla::real_t;
using namespace hyteg;

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

  P1Function< real_t > x("x", storage, minLevel, maxLevel);
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Edge>> emptyEdgeIds;
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Face>> emptyFaceIds;

  for (auto face : storage->getFaces())
  {
    for (size_t i = 0; i < v_perFace; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(face.second->getData(x.getFaceDataID())->getPointer(maxLevel)[i], 0.0);
    }
  }
  for (auto edge : storage->getEdges())
  {
    for (size_t i = 0; i < v_perEdge + edge.second->getNumHigherDimNeighbors() * nbr_v_perEdge; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(edge.second->getData(x.getEdgeDataID())->getPointer(maxLevel)[i], 0.0);
    }
  }
  for (auto vertex : storage->getVertices())
  {
    for (size_t i = 0; i < v_perVertex + vertex.second->getNumHigherDimNeighbors(); ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(vertex.second->getData(x.getVertexDataID())->getPointer(maxLevel)[i], 0.0);
    }
  }



  std::function<real_t(const Point3D &,const std::vector<real_t>&)> exact = [](const Point3D & xx, const std::vector<real_t>&) { return 2*xx[0] + xx[1]; };
  //std::function<real_t(const Point3D &)> exact = [](const Point3D & x) { return 13; };

  real_t value,xStepSize, yStepSize;
  for(auto faceIter : storage->getFaces())
  {
    auto face = faceIter.second;
    value = face->getCoordinates()[0][0] * 2 + face->getCoordinates()[0][0];
//    Edge *faceEdge0 = storage->getEdge(face->getEdgeID0().getID());
//    Edge *faceEdge1 = storage->getEdge(face->getEdgeID1().getID());
    xStepSize = walberla::real_c(face->getCoordinates()[1][0] - face->getCoordinates()[0][0]) / walberla::real_c((v_perEdge-1));
    yStepSize = walberla::real_c(face->getCoordinates()[2][1] - face->getCoordinates()[0][1]) / walberla::real_c((v_perEdge-1));

    vertexdof::macroface::interpolate< real_t >(maxLevel, *face, x.getFaceDataID(), emptyFaceIds, exact);
    for (uint_t i = 0; i < v_perEdge; ++i)
    {
      for (uint_t j = 0; j < v_perEdge - i; ++j)
      {
        uint_t idx = vertexdof::macroface::indexFromVertex( maxLevel, j, i,
                                                            stencilDirection::VERTEX_C );
        if (vertexdof::macroface::is_boundary(idx, v_perEdge))
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx], 0.0,
                                     "i: " << i << " j: " << j << " idx: " << idx << " value was " << value);
        } else
        {
          WALBERLA_CHECK_FLOAT_EQUAL(face->getData(x.getFaceDataID())->getPointer(maxLevel)[idx], value,
                                     "i: " << i << " j: " << j << " idx: " << idx << " value was " << value);
        }
        value += 2 * xStepSize;
      }
      value = yStepSize * (real_t) (i + 1);
    }
  }

  value = 0;
  for(auto edgeIter : storage->getEdges()){
    auto edge = edgeIter.second;
    hyteg::vertexdof::macroedge::interpolate< real_t >(maxLevel, *edge,x.getEdgeDataID(),emptyEdgeIds,exact);
    value = 2 * edge->getCoordinates()[0][0] + edge->getCoordinates()[0][1];
    xStepSize = edge->getDirection()[0] / walberla::real_c((v_perEdge-1));
    yStepSize = edge->getDirection()[1] / walberla::real_c((v_perEdge-1));
    value += 2*xStepSize;
    value += yStepSize;
    for(uint_t i = 1; i < v_perEdge-1; ++i)
    {
      WALBERLA_CHECK_FLOAT_EQUAL(edge->getData(x.getEdgeDataID())->getPointer(maxLevel)[i], value,
                                 "i: " << i << " edge: "<< *edge);
      value += 2*xStepSize;
      value += yStepSize;
    }
  }

  return EXIT_SUCCESS;
}
