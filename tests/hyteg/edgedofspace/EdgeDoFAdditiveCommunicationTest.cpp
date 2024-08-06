/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr.
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
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

static void
    testEdgeDoFAdditiveCommunication2D( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode,
                                        const std::string&                                                 meshFile )
{
   const uint_t level        = 3;
   const real_t testValue    = 1.0;
   const real_t someConstant = 6.345;

   auto                  meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupPrimitiveStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage );

   EdgeDoFFunction< real_t > x( "x", storage, level, level );
   x.setLocalCommunicationMode( localCommunicationMode );

   // 1. To avoid masking, interpolate any scalar != 0 everywhere
   x.interpolate( real_c( someConstant ), level, All );

   // 2. Set all entries (incl. ghost layers) on all cells to a test value.
   for ( const auto& it : storage->getFaces() )
   {
      auto facePtr = it.second->getData( x.getFaceDataID() )->getPointer( level );
      for ( const auto& idxIt : edgedof::macroface::Iterator( level, 0 ) )
      {
         for ( auto orientation : edgedof::faceLocalEdgeDoFOrientations )
         {
            const auto idx = edgedof::macroface::index( level, idxIt.x(), idxIt.y(), orientation );
            facePtr[idx]   = testValue;
         }
      }
   }

   // 3. Communicate additively. Each unknown on each primitive should now be equal to the number of neighbor faces times the test value.
   x.communicateAdditively< Face, Edge >( level, hyteg::DirichletBoundary, *storage );

   for ( const auto& it : storage->getEdges() )
   {
      auto edgePtr = it.second->getData( x.getEdgeDataID() )->getPointer( level );
      for ( const auto& idxIt : edgedof::macroedge::Iterator( level, 0 ) )
      {
         const auto idx = edgedof::macroedge::index( level, idxIt.x() );
         if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
         {
            //       this is always false but intel 2018 update 4 fails without
            if ( idx > 8 )
            {
               WALBERLA_LOG_DEVEL_VAR( idx );
            }
            WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], someConstant, "Face -> Edge additive comm failed (on Dirichlet boundary)" );
         }
         else
         {
            WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx],
                                        testValue * real_c( it.second->getNumNeighborFaces() ),
                                        "Face -> Edge additive comm failed (inner or Neumann)" );
         }
      }
   }
}

static void
    testEdgeDoFAdditiveCommunication3D( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode,
                                        const std::string&                                                 meshFile )
{
   const uint_t level        = 3;
   const real_t testValue    = 1.0;
   const real_t someConstant = 6.345;

   auto                  meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupPrimitiveStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage );

   EdgeDoFFunction< real_t > x( "x", storage, level, level );
   x.setLocalCommunicationMode( localCommunicationMode );

   // 1. To avoid masking, interpolate any scalar != 0 everywhere
   x.interpolate( real_c( someConstant ), level, All );

   // 2. Set all entries (incl. ghost layers) on all cells to a test value.
   for ( const auto& it : storage->getCells() )
   {
      auto cellPtr = it.second->getData( x.getCellDataID() )->getPointer( level );
      for ( const auto& idxIt : edgedof::macrocell::Iterator( level, 0 ) )
      {
         for ( auto orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
         {
            const auto idx = edgedof::macrocell::index( level, idxIt.x(), idxIt.y(), idxIt.z(), orientation );
            cellPtr[idx]   = testValue;
         }
      }
      for ( const auto& idxIt : edgedof::macrocell::IteratorXYZ( level, 0 ) )
      {
         const auto idx = edgedof::macrocell::index( level, idxIt.x(), idxIt.y(), idxIt.z(), edgedof::EdgeDoFOrientation::XYZ );
         cellPtr[idx]   = testValue;
      }
   }

   // 3. Communicate additively. Each unknown on each primitive should now be equal to the number of neighbor faces times the test value.
   x.communicateAdditively< Cell, Face >( level, hyteg::DirichletBoundary, *storage );

   for ( const auto& it : storage->getFaces() )
   {
      auto facePtr = it.second->getData( x.getFaceDataID() )->getPointer( level );
      for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
      {
         for ( const auto& idxIt : edgedof::macroface::Iterator( level, 0 ) )
         {
            if ( edgedof::macroface::isInnerEdgeDoF( level, idxIt, orientation ) )
            {
               const auto idx = edgedof::macroface::index( level, idxIt.x(), idxIt.y(), orientation );
               if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
               {
                  WALBERLA_CHECK_FLOAT_EQUAL(
                      facePtr[idx], someConstant, "Cell -> Face additive comm failed (on Dirichlet boundary)" );
               }
               else
               {
                  WALBERLA_CHECK_FLOAT_EQUAL( facePtr[idx],
                                              testValue * real_c( it.second->getNumNeighborCells() ),
                                              "Cell -> Face additive comm failed (inner or Neumann)" );
               }
            }
         }
      }
   }

   x.communicateAdditively< Cell, Edge >( level, hyteg::DirichletBoundary, *storage );

   for ( const auto& it : storage->getEdges() )
   {
      auto edgePtr = it.second->getData( x.getEdgeDataID() )->getPointer( level );
      for ( const auto& idxIt : edgedof::macroedge::Iterator( level, 0 ) )
      {
         const auto idx = edgedof::macroedge::index( level, idxIt.x() );
         if ( x.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
         {
            WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx], someConstant, "Cell -> Edge additive comm failed (on Dirichlet boundary)" );
         }
         else
         {
            WALBERLA_CHECK_FLOAT_EQUAL( edgePtr[idx],
                                        testValue * real_c( it.second->getNumNeighborCells() ),
                                        "Cell -> Edge additive comm failed (inner or Neumann)" );
         }
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "tri_1el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "tri_1el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "tri_2el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "tri_2el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "quad_8el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "quad_8el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "annulus_coarse.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication2D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "annulus_coarse.msh" ) );

   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ) );

   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "3D/pyramid_4el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "3D/pyramid_4el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                              hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ) );
   hyteg::testEdgeDoFAdditiveCommunication3D( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                              hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ) );

   return EXIT_SUCCESS;
}
