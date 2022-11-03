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
#include "hyteg/p1functionspace/VertexDoFPackInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "core/mpi/all.h"
#include "core/debug/all.h"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

int main (int argc, char ** argv )
{
  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);


  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  hyteg::P1Function< real_t > x("x", storage, minLevel, maxLevel);

  x.enumerate(maxLevel);
  hyteg::communication::syncFunctionBetweenPrimitives( x, maxLevel );

  uint_t numberOfChecks = 0;
  uint_t totalExpectedChecks = (3 * hyteg::levelinfo::num_microvertices_per_edge(maxLevel)
                                + 3 * ( hyteg::levelinfo::num_microvertices_per_edge(maxLevel) - 1))
                               * storage->getNumberOfLocalFaces();

  for (auto &faceIt : storage->getFaces()) {
    Face &face = *faceIt.second;
    //BubbleFace::printFunctionMemory<maxLevel>(face,x.getFaceDataID());

    real_t *faceData = face.getData(x.getFaceDataID())->getPointer(maxLevel);
    std::vector<PrimitiveID> nbrEdges;
    face.getNeighborEdges(nbrEdges);
    for(uint_t i = 0; i < nbrEdges.size(); ++i){
      Edge* edge = storage->getEdge(nbrEdges[0]);
      real_t* edgeData = edge->getData(x.getEdgeDataID())->getPointer(maxLevel);
      idx_t idxCounter = 0;
      uint_t faceIdOnEdge = edge->face_index(face.getID());
//////////////////// INNER VERTEX //////////////////////
      idxCounter = 0;
#if 0
      auto it = P1Face::indexIterator(face.edge_index(edge->getID()),
                                      face.edge_orientation[face.edge_index(edge->getID())],
                                      P1Face::VERTEX_INNER,
                                      maxLevel);
      for(; it != P1Face::indexIterator(); ++it){
#endif
      const indexing::FaceBoundaryDirection faceBorderDirection = indexing::getFaceBoundaryDirection(
          face.edge_index( edge->getID() ), face.getEdgeOrientation()[face.edge_index( edge->getID() )] );
      for ( const auto & it : vertexdof::macroface::BoundaryIterator( maxLevel, faceBorderDirection, 1 ) )
      {
        if(faceIdOnEdge == 0) {
          WALBERLA_CHECK_EQUAL(edgeData[vertexdof::macroedge::indexFromVertex( maxLevel, idxCounter, stencilDirection::VERTEX_SE )],
                               faceData[vertexdof::macroface::indexFromVertex(
          maxLevel, it.col(), it.row(), stencilDirection::VERTEX_C )]);
          numberOfChecks++;
        } else if(faceIdOnEdge == 1){
          WALBERLA_CHECK_EQUAL(edgeData[vertexdof::macroedge::indexFromVertex( maxLevel, idxCounter, stencilDirection::VERTEX_N )],
                               faceData[vertexdof::macroface::indexFromVertex(
          maxLevel, it.col(), it.row(), stencilDirection::VERTEX_C )]);
          numberOfChecks++;
        } else{
          WALBERLA_CHECK(false);
        }
        idxCounter++;
      }
//////////////////// VERTEX //////////////////////
      idxCounter = 0;
#if 0
      it = P1Face::indexIterator(face.edge_index(edge->getID()),
                                 face.edge_orientation[face.edge_index(edge->getID())],
                                 P1Face::VERTEX,
                                 maxLevel);
      for(; it != P1Face::indexIterator(); ++it){
#endif
      for ( const auto & it : vertexdof::macroface::BoundaryIterator( maxLevel, faceBorderDirection, 0 )  )
      {
        WALBERLA_CHECK_EQUAL(edgeData[vertexdof::macroedge::indexFromVertex( maxLevel, idxCounter, stencilDirection::VERTEX_C )],
                             faceData[vertexdof::macroface::indexFromVertex(
        maxLevel, it.col(), it.row(), stencilDirection::VERTEX_C )]);
        numberOfChecks++;
        idxCounter++;
      }
    }
  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks);


  numberOfChecks = 0;
  totalExpectedChecks = 0;
  for(auto &vertexIt : storage->getVertices())
  {
    //totalExpectedChecks++;
    totalExpectedChecks += vertexIt.second->getNumHigherDimNeighbors() * 2;
  }

  for (auto &edgeIt : storage->getEdges()) {

    Edge &edge = *edgeIt.second;
    //BubbleEdge::printFunctionMemory(edge,x.getEdgeDataID(),maxLevel);
    real_t *edgeData = edge.getData(x.getEdgeDataID())->getPointer(maxLevel);
    std::vector<PrimitiveID> nbrVertices;
    edge.getNeighborVertices(nbrVertices);
    for(uint_t i = 0; i < nbrVertices.size(); ++i)
    {
      Vertex* vertex = storage->getVertex(nbrVertices[i]);
      real_t* vertexData = vertex->getData(x.getVertexDataID())->getPointer(maxLevel);
      uint_t vertexIdOnEdge = edge.vertex_index(vertex->getID());
      uint_t vPerEdge = levelinfo::num_microvertices_per_edge(maxLevel);
      if(vertexIdOnEdge == 0){
        WALBERLA_CHECK_EQUAL(edgeData[vertexdof::macroedge::indexFromVertex( maxLevel, 0, stencilDirection::VERTEX_C )],
                             vertexData[0]);
        numberOfChecks++;

        WALBERLA_CHECK_EQUAL(edgeData[vertexdof::macroedge::indexFromVertex( maxLevel, 1, stencilDirection::VERTEX_C )],
                             vertexData[1 + vertex->edge_index(edge.getID())]);
        numberOfChecks++;

      } else if( vertexIdOnEdge == 1){
        WALBERLA_CHECK_EQUAL(edgeData[vertexdof::macroedge::indexFromVertex( maxLevel, idx_t( vPerEdge ) - 1, stencilDirection::VERTEX_C )],
                             vertexData[0]);
        numberOfChecks++;

        WALBERLA_CHECK_EQUAL(edgeData[vertexdof::macroedge::indexFromVertex( maxLevel, idx_t( vPerEdge ) - 2, stencilDirection::VERTEX_C )],
                             vertexData[1 + vertex->edge_index(edge.getID())]);
        numberOfChecks++;

      } else {
        WALBERLA_CHECK(false);
      }
    }
  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks);




  return 0;
}
