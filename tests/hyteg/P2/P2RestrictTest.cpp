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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

///this test check restrict on specific points on the grid
static void testP2Restrict() {
  const uint_t sourceLevel = 3;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");

  SetupPrimitiveStorage setupStorage(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  std::shared_ptr <PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  auto x = std::make_shared < P2Function < real_t > > ("x", storage, sourceLevel-1, sourceLevel, BoundaryCondition::create0123BC() );
  auto ident = P2Function < real_t >("ident", storage, sourceLevel-1, sourceLevel, BoundaryCondition::create0123BC() );
  std::function< real_t( const hyteg::Point3D& ) > one = []( const hyteg::Point3D&  ) {
     return 1;
  };
  ident.interpolate(one, sourceLevel, hyteg::All );

  x->enumerate( sourceLevel );
  x->add({1.0},{ident},sourceLevel, hyteg::All);
  hyteg::communication::syncP2FunctionBetweenPrimitives( *x, sourceLevel );

//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }
//
//    for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel>(*face, x->getEdgeDoFFunction()->getFaceDataID());
//  }
//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

  P2toP2QuadraticRestriction restrictionOperator;
  restrictionOperator.restrict( *x, sourceLevel, All );



//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel - 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }
//
//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel-1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

  ///calculate expected entry for idx = 1 on edge 0
  real_t expected = 5;
  expected += (1./8.) * (-46 + 3 * 47 + 3 * 48 - 49);
  expected += (1./8.) * (-98 + 3 * 99 + 3 *127 - 128);
  expected += (1./8.) * (-70 - 72);
  expected += (1./8.) * (-105 - 133 - 106 - 134);

  WALBERLA_CHECK_FLOAT_EQUAL(storage->getEdge(PrimitiveID(3))->getData(x->getVertexDoFFunction().getEdgeDataID())->getPointer(sourceLevel - 1)[1], expected);


  ///calculate expected entry for idx = 2 on edge 1
  expected = 21;
  expected += (1./8.) * (- 64 + 3 * 65 + 3 * 66 - 67);
  expected += (1./8.) * (- 143 + 3 * 147 + 3 * 91 - 94);
  expected += (1./8.) * (- 115 - 122);
  expected += (1./8.) * (- 142 - 146 - 86 - 90);

  WALBERLA_CHECK_FLOAT_EQUAL(storage->getEdge(PrimitiveID(5))->getData(x->getVertexDoFFunction().getEdgeDataID())->getPointer(sourceLevel - 1)[2], expected);

  ///calculate expected entry for idx = 2 on edge 1
  expected = 16;
  expected += (1./8.) * (- 61 + 3 * 60 + 3 * 59 - 58);
  expected += (1./8.) * (- 120 + 3 * 123 + 3 * 95 - 97);
  expected += (1./8.) * (- 148 - 153);
  expected += (1./8.) * (- 124 - 121 - 96 - 93);

  WALBERLA_CHECK_FLOAT_EQUAL(storage->getEdge(PrimitiveID(4))->getData(x->getVertexDoFFunction().getEdgeDataID())->getPointer(sourceLevel - 1)[3], expected);

  ///calculate expected entry for vertex dof on face at 4,2
  expected = 34;
  /// horizontal edges
  expected += (1./8.) * (- 85 - 87);
  expected += (1./8.) * (- 79 + 3 * 80 + 3 * 81 - 82);
  expected += (1./8.) * (- 73 - 75);
  /// vertical edges
  expected += (1./8.) * (- 146 - 147);
  expected += (1./8.) * (3 * 142 - 143);
  expected += (1./8.) * (-135 + 3 * 136);
  expected += (1./8.) * (-129 - 130);
  /// diagonal edges
  expected += (1./8.) * (- 118 - 119);
  expected += (1./8.) * (- 113 + 3 * 114);
  expected += (1./8.) * ( 3 * 109- 110);
  expected += (1./8.) * (- 102 - 103);
  WALBERLA_CHECK_FLOAT_EQUAL(storage->
    getFace(PrimitiveID(6))->
    getData(x->getVertexDoFFunction().getFaceDataID())->
    getPointer(sourceLevel - 1)[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel - 1,2,1,stencilDirection::VERTEX_C)], expected);


}

/// sets all values on the source level to a specific value and checks after one restrict step
static void testP2Restrict2() {
  const uint_t sourceLevel = 4;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = std::make_shared<P2Function<real_t> >("x", storage, sourceLevel - 2, sourceLevel, BoundaryCondition::create0123BC() );
  typedef stencilDirection sD;
  std::function<real_t(const hyteg::Point3D &)> values = [](const hyteg::Point3D &) { return 13; };

  x->interpolate(values, sourceLevel);

  P2toP2QuadraticRestriction restrictionOperator;
  restrictionOperator.restrict( *x, sourceLevel, All );

//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel - 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }
//
//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel - 1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }
//
//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel - 1>(*face, x->getEdgeDoFFunction()->getFaceDataID());
//  }
//
//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::edgedof::macroedge::printFunctionMemory<real_t, sourceLevel - 1>(*edge, x->getEdgeDoFFunction()->getEdgeDataID());
//  }

  real_t* edgeDoFCoarseData = storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(sourceLevel - 1);
  real_t* vertexDoFCoarseData = storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction().getFaceDataID())->getPointer(sourceLevel - 1);

  for( const auto & it : hyteg::vertexdof::macroface::Iterator( sourceLevel - 1, 1)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel - 1,it.col(), it.row(), sD::VERTEX_C)],
      13.,
      it.col() << " " << it.row());
  }

  for( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel - 1, 0)) {
    if(it.row() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel - 1,it.col(), it.row(), sD::EDGE_HO_E)],
        65.,
        it.col() << " " << it.row());
    }
    if(it.col() + it.row() != idx_t( hyteg::levelinfo::num_microedges_per_edge( sourceLevel - 1 ) - 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel - 1,it.col(), it.row(), sD::EDGE_DI_NE)],
        65.,
        it.col() << " " << it.row());
    }
    if(it.col() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel - 1,it.col(), it.row(), sD::EDGE_VE_N)],
        65.,
        it.col() << " " << it.row());
    }
  }

  restrictionOperator.restrict( *x, sourceLevel - 1, All );

  edgeDoFCoarseData = storage->getFace(PrimitiveID(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(sourceLevel - 2);
  vertexDoFCoarseData = storage->getFace(PrimitiveID(6))->getData(x->getVertexDoFFunction().getFaceDataID())->getPointer(sourceLevel - 2);

  for( const auto & it : hyteg::vertexdof::macroface::Iterator( sourceLevel - 2, 1)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFCoarseData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel - 2,it.col(), it.row(), sD::VERTEX_C)],
      13.,
      it.col() << " " << it.row());
  }

  for( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel - 2, 0)) {
    if(it.row() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel - 2,it.col(), it.row(), sD::EDGE_HO_E)],
        273.,
        it.col() << " " << it.row());
    }
    if(it.col() + it.row() != idx_t( hyteg::levelinfo::num_microedges_per_edge( sourceLevel - 2 ) - 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel - 2,it.col(), it.row(), sD::EDGE_DI_NE)],
        273.,
        it.col() << " " << it.row());
    }
    if(it.col() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFCoarseData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel - 2,it.col(), it.row(), sD::EDGE_VE_N)],
        273.,
        it.col() << " " << it.row());
    }
  }

}


static void testP2InterpolateAndRestrict() {

  const uint_t sourceLevel = 5;
  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = P2Function<real_t> ("x", storage, sourceLevel - 3, sourceLevel);
  auto y = P2Function<real_t> ("y", storage, sourceLevel - 3, sourceLevel);
  auto error = P2Function<real_t> ("x", storage, sourceLevel - 3, sourceLevel);

  std::function<real_t(const hyteg::Point3D&)> exact = [](const hyteg::Point3D& xx) { return 2. * xx[0] * xx[0] + 3. * xx[0] + 13. + 4. * xx[1] + 5. *  xx[1] * xx[1]; };
  x.interpolate(exact,sourceLevel    , hyteg::All);
  y.interpolate(exact,sourceLevel - 3, hyteg::All);

  x.restrictInjection(sourceLevel    , hyteg::All);
  x.restrictInjection(sourceLevel - 1, hyteg::All);
  x.restrictInjection(sourceLevel - 2, hyteg::All);

  error.assign({1.0, -1.0}, {x, y}, sourceLevel - 3, hyteg::All);

  WALBERLA_CHECK_FLOAT_EQUAL(error.dotGlobal(error,sourceLevel - 3, hyteg::All),0.);
}

}/// namespace hyteg

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hyteg::testP2Restrict();
  hyteg::testP2Restrict2();
  hyteg::testP2InterpolateAndRestrict();

  return EXIT_SUCCESS;
}
