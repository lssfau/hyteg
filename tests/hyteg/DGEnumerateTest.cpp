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
#include "core/debug/all.h"
#include "core/mpi/all.h"

#include "hyteg/facedofspace/FaceDoFFunction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using namespace hyteg;

using walberla::real_t;

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

  hyteg::FaceDoFFunction< uint_t > x("x", storage, minLevel, maxLevel);

  uint_t check = 1;
  uint_t sum = 0;
  x.enumerate(maxLevel,1);

  uint_t totalDoFs = ( hyteg::levelinfo::num_microfaces_per_face(maxLevel) * storage->getNumberOfLocalFaces());
  uint_t expectedSum = (totalDoFs * (totalDoFs + 1))/2;

  for (auto &vertexIt : storage->getVertices()){
    Vertex &vertex = *vertexIt.second;
    uint_t *vertexData = vertex.getData(x.getVertexDataID())->getPointer(maxLevel);
    for(uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i){
      WALBERLA_CHECK_EQUAL(vertexData[i*2],check);
      sum += check;
      check++;
    }
  }



  for (auto &edgeIt : storage->getEdges()) {
    Edge &edge = *edgeIt.second;
    uint_t *edgeData = edge.getData(x.getEdgeDataID())->getPointer(maxLevel);
    uint_t FaceDoFonFace = hyteg::levelinfo::num_microvertices_per_edge(maxLevel) * 2 - 3 ;
    //this only works with the linear default layout; can be changed to use index function
    for(uint_t i = 0; i < edge.getNumHigherDimNeighbors(); ++i){
      uint_t start = FaceDoFonFace * i;
      for(uint_t j = 2; j < hyteg::levelinfo::num_microvertices_per_edge(maxLevel) * 2 - 5; j += 2){
        WALBERLA_CHECK_EQUAL(edgeData[start + j],check);
        sum += check;
        check++;
      }
    }
  }

  for ( auto &faceIt : storage->getFaces() ) {
    Face &face = *faceIt.second;
    uint_t *faceData = face.getData(x.getFaceDataID())->getPointer(maxLevel);
    uint_t rowsize = levelinfo::num_microvertices_per_edge(maxLevel) - 2;
    uint_t inner_rowsize = rowsize;
    for(uint_t i = 1; i < (rowsize -1 ); ++i){
      for(uint_t j = 1; j < ( inner_rowsize - 1 ); ++j){
        WALBERLA_CHECK_EQUAL(
          faceData[ facedof::macroface::indexFaceFromVertex( maxLevel, i, j, stencilDirection::CELL_GRAY_NE ) ],
          check);
        sum += check;
        ++check;
      }
      --inner_rowsize;
    }
    inner_rowsize = rowsize;
    for(uint_t i = 0; i < rowsize ; ++i){
      for(uint_t j = 0; j <  inner_rowsize ; ++j){
        WALBERLA_CHECK_EQUAL(
          faceData[ facedof::macroface::indexFaceFromGrayFace( maxLevel, i, j, stencilDirection::CELL_BLUE_E ) ],
          check);
        sum += check;
        ++check;
      }
      --inner_rowsize;
    }
  }
  --check;
  WALBERLA_CHECK_EQUAL(check,totalDoFs);
  WALBERLA_CHECK_EQUAL(sum,expectedSum);

  return 0;
}
