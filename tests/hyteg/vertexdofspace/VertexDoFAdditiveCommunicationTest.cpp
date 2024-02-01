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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

namespace hyteg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;


static void testVertexDoFAdditiveCommunication2D( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode,
                                                  const std::string & meshFile )
{
  const uint_t level = 3;
  const real_t testValue = 1.0;
  const real_t someConstant = 6.345;


  auto meshInfo = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupPrimitiveStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  auto storage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage );

  vertexdof::VertexDoFFunction< real_t > x( "x", storage, level, level );
  x.setLocalCommunicationMode( localCommunicationMode );

  // 1. To avoid masking, interpolate any scalar != 0 everywhere
  x.interpolate( real_c( someConstant ), level, All );

  // 2. Set all entries (incl. ghost layers) on all cells to a test value.
  for ( const auto & it : storage->getFaces() )
  {
    auto facePtr = it.second->getData( x.getFaceDataID() )->getPointer( level );
    for ( const auto & idxIt : vertexdof::macroface::Iterator( level, 0 ) )
    {
      const auto idx = vertexdof::macroface::indexFromVertex( level, idxIt.x(), idxIt.y(), stencilDirection::VERTEX_C );
      facePtr[ idx ] = testValue;
    }
  }

  // 3. Communicate additively. Each unknown on each primitive should now be equal to the number of neighbor faces times the test value.
  x.communicateAdditively< Face, Edge >( level, hyteg::DirichletBoundary, *storage );
  x.communicateAdditively< Face, Vertex >( level, hyteg::DirichletBoundary, *storage );

  for ( const auto & it : storage->getEdges() )
  {
    auto edgePtr = it.second->getData( x.getEdgeDataID() )->getPointer( level );
    for ( const auto & idxIt : vertexdof::macroedge::Iterator( level, 1 ) )
    {
      const auto idx = vertexdof::macroedge::indexFromVertex( level, idxIt.x(), stencilDirection::VERTEX_C );
      if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], someConstant, "Face -> Edge additive comm failed (on Dirichlet boundary)" );
      }
      else
      {
        WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], testValue * real_c( it.second->getNumNeighborFaces() ), "Face -> Edge additive comm failed (inner or Neumann)" );
      }
    }
  }

  for ( const auto & it : storage->getVertices() )
  {
    auto vertexPtr = it.second->getData( x.getVertexDataID() )->getPointer( level );
    if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( vertexPtr[0], someConstant, "Face -> Vertex additive comm failed (on Dirichlet boundary)" );
    }
    else
    {
      WALBERLA_CHECK_FLOAT_EQUAL( vertexPtr[0], testValue * real_c( it.second->getNumNeighborFaces() ), "Face -> Vertex additive comm failed (inner or Neumann)" );
    }
  }

}



static void testVertexDoFAdditiveCommunication3D( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode,
                                                  const std::string & meshFile )
{
  const uint_t level = 3;
  const real_t testValue = 1.0;
  const real_t someConstant = 6.345;


  auto meshInfo = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupPrimitiveStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  auto storage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage );

  vertexdof::VertexDoFFunction< real_t > x( "x", storage, level, level );
  x.setLocalCommunicationMode( localCommunicationMode );

  // 1. To avoid masking, interpolate any scalar != 0 everywhere
  x.interpolate( real_c( someConstant ), level, All );

  // 2. Set all entries (incl. ghost layers) on all cells to a test value.
  for ( const auto & it : storage->getCells() )
  {
    auto cellPtr = it.second->getData( x.getCellDataID() )->getPointer( level );
    for ( const auto & idxIt : vertexdof::macrocell::Iterator( level, 0 ) )
    {
      const auto idx = vertexdof::macrocell::indexFromVertex( level, idxIt.x(), idxIt.y(), idxIt.z(), stencilDirection::VERTEX_C );
      cellPtr[ idx ] = testValue;
    }
  }

  // 3. Communicate additively. Each unknown on each primitive should now be equal to the number of neighbor cells times the test value.
  x.communicateAdditively< Cell, Face >( level, hyteg::DirichletBoundary, *storage );
  x.communicateAdditively< Cell, Edge >( level, hyteg::DirichletBoundary, *storage );
  x.communicateAdditively< Cell, Vertex >( level, hyteg::DirichletBoundary, *storage );

  for ( const auto & it : storage->getFaces() )
  {
    auto facePtr = it.second->getData( x.getFaceDataID() )->getPointer( level );

    for ( const auto & idxIt : vertexdof::macroface::Iterator( level, 1 ) )
    {
      const auto idx = vertexdof::macroface::indexFromVertex( level, idxIt.x(), idxIt.y(), stencilDirection::VERTEX_C );
      if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( facePtr[idx], someConstant, "Cell -> Face additive comm failed (on Dirichlet boundary)" );
      }
      else
      {
        WALBERLA_CHECK_FLOAT_EQUAL( facePtr[idx], testValue * real_c( it.second->getNumNeighborCells() ), "Cell -> Face additive comm failed (inner or Neumann)" );
      }
    }
  }

  for ( const auto & it : storage->getEdges() )
  {
    auto edgePtr = it.second->getData( x.getEdgeDataID() )->getPointer( level );
    for ( const auto & idxIt : vertexdof::macroedge::Iterator( level, 1 ) )
    {
      const auto idx = vertexdof::macroedge::indexFromVertex( level, idxIt.x(), stencilDirection::VERTEX_C );
      if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], someConstant, "Cell -> Edge additive comm failed (on Dirichlet boundary)" );
      }
      else
      {
        WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], testValue * real_c( it.second->getNumNeighborCells() ), "Cell -> Edge additive comm failed (inner or Neumann)" );
      }
    }
  }

  for ( const auto & it : storage->getVertices() )
  {
    auto vertexPtr = it.second->getData( x.getVertexDataID() )->getPointer( level );
    if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( vertexPtr[0], someConstant, "Cell -> Vertex additive comm failed (on Dirichlet boundary)" );
    }
    else
    {
      WALBERLA_CHECK_FLOAT_EQUAL( vertexPtr[0], testValue * real_c( it.second->getNumNeighborCells() ), "Cell -> Vertex additive comm failed (inner or Neumann)" );
    }
  }

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/tri_1el.msh" );
  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/tri_1el.msh" );
  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/tri_2el.msh" );
  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/tri_2el.msh" );
  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/quad_8el.msh" );
  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/quad_8el.msh" );
  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/annulus_coarse.msh" );
  hyteg::testVertexDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,       "../../data/meshes/annulus_coarse.msh" );

  hyteg::testVertexDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/3D/tet_1el.msh" );
  hyteg::testVertexDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT      , "../../data/meshes/3D/tet_1el.msh" );
  hyteg::testVertexDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/3D/pyramid_4el.msh" );
  hyteg::testVertexDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT      , "../../data/meshes/3D/pyramid_4el.msh" );
  hyteg::testVertexDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI, "../../data/meshes/3D/regular_octahedron_8el.msh" );
  hyteg::testVertexDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT      , "../../data/meshes/3D/regular_octahedron_8el.msh" );

  return EXIT_SUCCESS;
}
