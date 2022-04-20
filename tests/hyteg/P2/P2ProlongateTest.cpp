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

#include "hyteg/StencilDirections.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

/// this test checks if a constant value is prolongated correctly to all points on the finer level
static void testP2Prolongate() {
  const uint_t sourceLevel = 2;

  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = std::make_shared<P2Function<real_t> >("x", storage, sourceLevel, sourceLevel + 2, BoundaryCondition::create0123BC());
  typedef stencilDirection sD;
  std::function<real_t(const hyteg::Point3D &)> values = [](const hyteg::Point3D &) { return 13; };
  x->interpolate(values, sourceLevel, hyteg::All);

  P2toP2QuadraticProlongation prolongationOperator;
  prolongationOperator.prolongate( *x, sourceLevel, All );

  real_t* edgeDoFFineData = storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(sourceLevel + 1);
  real_t* vertexDoFFineData = storage->getFace(PrimitiveID::create(6))->getData(x->getVertexDoFFunction().getFaceDataID())->getPointer(sourceLevel + 1);

  for( const auto & it : hyteg::vertexdof::macroface::Iterator( sourceLevel + 1, 1)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel + 1,it.col(), it.row(), sD::VERTEX_C)],
      13.,
      it.col() << " " << it.row());
  }

  for( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel + 1, 0)) {
    if(it.row() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,it.col(), it.row(), sD::EDGE_HO_E)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( sourceLevel + 1 ) - 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,it.col(), it.row(), sD::EDGE_DI_NE)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,it.col(), it.row(), sD::EDGE_VE_N)],
        13.,
        it.col() << " " << it.row());
    }
  }


  for (auto &edgeIT : storage->getEdges()) {
    auto edge = edgeIT.second;
    edgeDoFFineData = edge->getData(x->getEdgeDoFFunction().getEdgeDataID())->getPointer(
      sourceLevel + 1);
    vertexDoFFineData = edge->getData(x->getVertexDoFFunction().getEdgeDataID())->getPointer(
      sourceLevel + 1);

    for (const auto &it : hyteg::vertexdof::macroedge::Iterator(sourceLevel + 1, 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        vertexDoFFineData[hyteg::vertexdof::macroedge::indexFromVertex(sourceLevel + 1,it.col(), sD::VERTEX_C)],
        13.,
        it.col() << " " << it.row());
    }

    for (const auto &it : hyteg::edgedof::macroedge::Iterator(sourceLevel + 1, 0)) {
        WALBERLA_CHECK_FLOAT_EQUAL(
          edgeDoFFineData[hyteg::edgedof::macroedge::indexFromVertex(sourceLevel + 1,it.col(), sD::EDGE_HO_E)],
          13.,
          it.col() << " " << it.row());
    }
  }

  prolongationOperator.prolongate( *x, sourceLevel + 1, All );

  edgeDoFFineData = storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(sourceLevel + 2);
  vertexDoFFineData = storage->getFace(PrimitiveID::create(6))->getData(x->getVertexDoFFunction().getFaceDataID())->getPointer(sourceLevel + 2);

  for( const auto & it : hyteg::vertexdof::macroface::Iterator( sourceLevel + 2, 1)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel + 2,it.col(), it.row(), sD::VERTEX_C)],
      13.,
      it.col() << " " << it.row());
  }

  for( const auto & it : hyteg::edgedof::macroface::Iterator( sourceLevel + 2, 0)) {
    if(it.row() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 2,it.col(), it.row(), sD::EDGE_HO_E)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( sourceLevel + 2 ) - 1)) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 2,it.col(), it.row(), sD::EDGE_DI_NE)],
        13.,
        it.col() << " " << it.row());
    }
    if(it.col() != 0) {
      WALBERLA_CHECK_FLOAT_EQUAL(
        edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 2,it.col(), it.row(), sD::EDGE_VE_N)],
        13.,
        it.col() << " " << it.row());
    }
  }


//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getEdgeDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::edgedof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getEdgeDoFFunction()->getEdgeDataID());
//  }

}


/// this test writes specific values at certain points and
/// checks wether these values are propagated correctly after one prolongation step
static void testP2Prolongate2() {

  const uint_t sourceLevel = 2;
  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/tri_1el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = std::make_shared<P2Function<real_t> >("x", storage, sourceLevel, sourceLevel + 1, BoundaryCondition::create0123BC());
  typedef stencilDirection sD;

  P2toP2QuadraticProlongation prolongationOperator;

  /// this should not be necessary but just to be save
  std::function<real_t(const hyteg::Point3D &)> zeros = [](const hyteg::Point3D &) { return 0; };
  x->interpolate(zeros, sourceLevel);

  real_t* edgeDoFFineData = storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(sourceLevel + 1);
  real_t* vertexDoFFineData = storage->getFace(PrimitiveID::create(6))->getData(x->getVertexDoFFunction().getFaceDataID())->getPointer(sourceLevel + 1);

  /// all possible vertical edge Dof locations that need to be updated by the face
  std::vector<std::pair<idx_t, idx_t > > vertical   = { {1,0},{2,0},{3,0},{1,1},{2,1},{1,2} };
  /// all possible horizontal edge Dof locations that need to be updated by the face
  std::vector<std::pair<idx_t, idx_t > > horizontal = { {0,1},{1,1},{2,1},{0,2},{1,2},{0,3} };
  /// all possible diagonal edge Dof locations that need to be updated by the face
  std::vector<std::pair<idx_t, idx_t > > diagonal   = { {0,0},{1,0},{2,0},{0,1},{1,1},{0,2} };
  /// all possible vertex dof locations that need to be updated by the face
  std::vector<std::pair<idx_t, idx_t > > vertex     = { {1,1},{2,1},{1,2} };

///////////////////////////
/// CHECH VERTICAL EDGE ///
///////////////////////////
  for(auto p : vertical) {

    storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::EDGE_VE_N)] = 16.0;

    prolongationOperator.prolongate( *x, sourceLevel, All );

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_VE_N)], 12.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_VE_S)], 12.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_HO_E)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_HO_W)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_DI_NW)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_DI_SE)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_VE_SE)], 4.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::EDGE_VE_NW)], 4.);

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2, p.second *2 + 1, sD::VERTEX_C)], 16.);

    storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::EDGE_VE_N)] = 0.0;

  }
/////////////////////////////
/// CHECK HORIZONTAL EDGE ///
/////////////////////////////
  for(auto p : horizontal) {

    storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::EDGE_HO_E)] = 16.0;

    prolongationOperator.prolongate( *x, sourceLevel, All );

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_VE_N)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_VE_S)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_HO_E)], 12.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_HO_W)], 12.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_DI_NW)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_DI_SE)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_HO_SE)], 4.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::EDGE_HO_NW)], 4.);

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2, sD::VERTEX_C)], 16.);

    storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::EDGE_HO_E)] = 0.0;

  }
///////////////////////////
/// CHECK DIAGONAL EDGE ///
///////////////////////////
  for(auto p : diagonal) {

    storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::EDGE_DI_NE)] = 16.0;

    prolongationOperator.prolongate( *x, sourceLevel, All );

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_VE_N)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_VE_S)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_HO_E)], 8.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_HO_W)], 8.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_NW)], 12.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_SE)], 12.);

    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_SW)], 4.);
    WALBERLA_CHECK_FLOAT_EQUAL(edgeDoFFineData[hyteg::edgedof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::EDGE_DI_NE)], 4.);

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 + 1, p.second * 2 + 1, sD::VERTEX_C)], 16.);

    storage->getFace(PrimitiveID::create(6))->getData(x->getEdgeDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::edgedof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::EDGE_DI_NE)] = 0.0;

  }
////////////////////
/// CHECK VERTEX ///
////////////////////
  for( auto p : vertex ){
    storage->getFace(PrimitiveID::create(6))->getData(x->getVertexDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::VERTEX_C)] = 16.0;

    prolongationOperator.prolongate( *x, sourceLevel, All );

    WALBERLA_CHECK_FLOAT_EQUAL(vertexDoFFineData[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel + 1,p.first * 2 , p.second * 2 , sD::VERTEX_C)],
                               16.,
                               p.first << " " << p.second);

    storage->getFace(PrimitiveID::create(6))->getData(x->getVertexDoFFunction().getFaceDataID())->getPointer(
      sourceLevel)[hyteg::vertexdof::macroface::indexFromVertex(sourceLevel,p.first, p.second, stencilDirection::VERTEX_C)] = 0.0;

  }
///////////////////////////



//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::vertexdof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getVertexDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::vertexdof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getVertexDoFFunction()->getEdgeDataID());
//  }

//  for (auto &faceIT : storage->getFaces()) {
//    auto face = faceIT.second;
//    hyteg::edgedof::macroface::printFunctionMemory<real_t, sourceLevel + 1>(*face, x->getEdgeDoFFunction()->getFaceDataID());
//  }

//  for (auto &edgeIT : storage->getEdges()) {
//    auto edge = edgeIT.second;
//    hyteg::edgedof::macroedge::printFunctionMemory<real_t, sourceLevel + 1>(*edge, x->getEdgeDoFFunction()->getEdgeDataID());
//  }

}

static void testP2InterpolateAndProlongate() {

  const uint_t sourceLevel = 2;
  const uint_t targetLevel = sourceLevel + 3;
  MeshInfo mesh = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
  std::shared_ptr<SetupPrimitiveStorage> setupStorage =
    std::make_shared<SetupPrimitiveStorage>(mesh, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(*setupStorage);
  auto x = P2Function<real_t> ("x", storage, sourceLevel, targetLevel, BoundaryCondition::create0123BC());
  auto y = P2Function<real_t> ("y", storage, sourceLevel, targetLevel, BoundaryCondition::create0123BC());
  auto error = P2Function<real_t> ("x", storage, sourceLevel, targetLevel, BoundaryCondition::create0123BC());

  std::function<real_t(const hyteg::Point3D&)> exact = [](const hyteg::Point3D& xx) { return 2. * xx[0] * xx[0] + 3. * xx[0] + 13. + 4. * xx[1] + 5. *  xx[1] * xx[1]; };
  x.interpolate(exact,sourceLevel    , hyteg::All);
  y.interpolate(exact,targetLevel, hyteg::All);

  P2toP2QuadraticProlongation prolongationOperator;

  prolongationOperator.prolongate( x, sourceLevel    , All );
  prolongationOperator.prolongate( x, sourceLevel + 1, All );
  prolongationOperator.prolongate( x, sourceLevel + 2, All );

  error.assign({1.0, -1.0}, {x, y}, targetLevel, hyteg::All);

  WALBERLA_CHECK_FLOAT_EQUAL(error.dotGlobal(error,targetLevel, hyteg::All),0.);
}

}/// namespace hyteg

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();
  hyteg::testP2Prolongate();
  hyteg::testP2Prolongate2();
  hyteg::testP2InterpolateAndProlongate();

  return EXIT_SUCCESS;
}
