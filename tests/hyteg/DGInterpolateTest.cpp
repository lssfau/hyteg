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
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/facedofspace_old/FaceDoFFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::real_c;
using namespace hyteg;

int main(int argc, char **argv) {
  walberla::debug::enterTestMode();
  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/quad_2el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);


  const size_t minLevel = 2;
  const size_t maxLevel = 5;

  size_t v_perEdge = levelinfo::num_microvertices_per_edge(maxLevel);

  std::function<real_t(const hyteg::Point3D&,const std::vector<real_t>&)> exact =
    [&](const hyteg::Point3D& xx,const std::vector<real_t>&) { return 2*xx[0] + xx[1]; };


  std::shared_ptr<Edge> edgeWithTwoFaces;
  for(auto edgeIt : storage->getEdges()){
    if(edgeIt.second.get()->getNumNeighborFaces() == 2){
      edgeWithTwoFaces = edgeIt.second;
    }
  }
  FaceDoFFunction_old< real_t > x("x", storage, minLevel, maxLevel);

  facedof::macroedge::interpolate< real_t >(maxLevel, *edgeWithTwoFaces,x.getEdgeDataID(),{},exact,storage);
  //facedof::macroedge::interpolate< real_t >(maxLevel,*edgeWithTwoFaces,x.getEdgeDataID(),{},exact,storage);


  Point3D edgeZeroPoint = edgeWithTwoFaces->getCoordinates()[0];
  Point3D edgeDirectionDx = edgeWithTwoFaces->getDirection()/(real_c(v_perEdge) -1.);
  auto face0 = storage->getFace(edgeWithTwoFaces->neighborFaces()[0]);
  uint_t vertexOnFace0 = face0->vertex_index(face0->get_vertex_opposite_to_edge(edgeWithTwoFaces->getID()));
  Point3D face0Dx = (face0->getCoordinates()[vertexOnFace0] - edgeZeroPoint) / (real_c(v_perEdge) - 1.);
  Point3D xPoint = edgeZeroPoint + 1.0/3.0 * (face0Dx + edgeDirectionDx) + edgeDirectionDx;
  for(uint_t i = 1; i < v_perEdge - 2; ++i){
    WALBERLA_CHECK_FLOAT_EQUAL(
      edgeWithTwoFaces->getData(x.getEdgeDataID())->getPointer(maxLevel)[facedof::macroedge::indexFaceFromVertex( maxLevel, i, stencilDirection::CELL_GRAY_SE )],
      exact(xPoint,{}), i )
    xPoint += edgeDirectionDx;
  }

  auto face1 = storage->getFace(edgeWithTwoFaces->neighborFaces()[1]);
  uint_t vertexOnFace1 = face1->vertex_index(face1->get_vertex_opposite_to_edge(edgeWithTwoFaces->getID()));
  Point3D face1Dx = (face1->getCoordinates()[vertexOnFace1] - edgeZeroPoint) / (real_c(v_perEdge) - 1.);
  xPoint = edgeZeroPoint + 1.0/3.0 * (face1Dx + edgeDirectionDx) + edgeDirectionDx;
  for(uint_t i = 1; i < v_perEdge - 2; ++i){
    WALBERLA_CHECK_FLOAT_EQUAL(
      edgeWithTwoFaces->getData(x.getEdgeDataID())->getPointer(maxLevel)[facedof::macroedge::indexFaceFromVertex( maxLevel, i, stencilDirection::CELL_GRAY_NE )],
      exact(xPoint,{}), i )
    xPoint += edgeDirectionDx;
  }

  //BubbleEdge::printFunctionMemory(*edgeWithTwoFaces,x.getEdgeDataID(),maxLevel);
}
