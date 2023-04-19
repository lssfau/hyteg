/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

static void
    testN1E1AdditiveCommunication3D( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode,
                                     const std::string&                                                 meshFile )
{
   const uint_t  level        = 3;
   const Point3D testValue    = { 1.0, 2.0, 3.0 };
   const real_t  someConstant = 6.345;

   auto                  meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupPrimitiveStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupPrimitiveStorage );

   n1e1::N1E1VectorFunction< real_t > ref( "ref", storage, level, level );
   n1e1::N1E1VectorFunction< real_t > test( "test", storage, level, level );
   test.setLocalCommunicationMode( localCommunicationMode );

   // 1. Interpolate a constant everywhere.
   ref.interpolate( testValue, level, All );
   test.interpolate( testValue, level, All );

   // 2. Set DoFs on lower dim primitives to some constant.
   test.getDoFs()->interpolateByPrimitiveType< Face >( someConstant, level, All );
   test.getDoFs()->interpolateByPrimitiveType< Edge >( someConstant, level, All );

   // 3. Communicate additively.
   // Each unknown on each primitive should now be equal to the number of neighbor faces times the reference value.
   test.communicateAdditively< Cell, Face >( level, hyteg::DirichletBoundary, *storage );

   for ( const auto& it : storage->getFaces() )
   {
      auto refFacePtr  = it.second->getData( ref.getDoFs()->getFaceDataID() )->getPointer( level );
      auto testFacePtr = it.second->getData( test.getDoFs()->getFaceDataID() )->getPointer( level );
      for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
      {
         for ( const auto& idxIt : edgedof::macroface::Iterator( level, 0 ) )
         {
            if ( edgedof::macroface::isInnerEdgeDoF( level, idxIt, orientation ) )
            {
               const auto idx = edgedof::macroface::index( level, idxIt.x(), idxIt.y(), orientation );
               if ( test.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) ==
                    DoFType::DirichletBoundary )
               {
                  WALBERLA_CHECK_FLOAT_EQUAL(
                      testFacePtr[idx], someConstant, "Cell -> Face additive comm failed (on Dirichlet boundary)" );
               }
               else
               {
                  WALBERLA_CHECK_FLOAT_EQUAL( testFacePtr[idx],
                                              refFacePtr[idx] * real_c( it.second->getNumNeighborCells() ),
                                              "Cell -> Face additive comm failed (inner or Neumann)" );
               }
            }
         }
      }
   }

   test.communicateAdditively< Cell, Edge >( level, hyteg::DirichletBoundary, *storage );

   for ( const auto& it : storage->getEdges() )
   {
      auto refEdgePtr  = it.second->getData( ref.getDoFs()->getEdgeDataID() )->getPointer( level );
      auto testEdgePtr = it.second->getData( test.getDoFs()->getEdgeDataID() )->getPointer( level );
      for ( const auto& idxIt : edgedof::macroedge::Iterator( level, 0 ) )
      {
         const auto idx = edgedof::macroedge::index( level, idxIt.x() );
         if ( test.getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() ) == DoFType::DirichletBoundary )
         {
            WALBERLA_CHECK_FLOAT_EQUAL(
                testEdgePtr[idx], someConstant, "Cell -> Edge additive comm failed (on Dirichlet boundary)" );
         }
         else
         {
            WALBERLA_CHECK_FLOAT_EQUAL( testEdgePtr[idx],
                                        refEdgePtr[idx] * real_c( it.second->getNumNeighborCells() ),
                                        "Cell -> Edge additive comm failed (inner or Neumann)" );
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( meshFile << ": passed" )
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                    "../../data/meshes/3D/tet_1el.msh" );
   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                    "../../data/meshes/3D/tet_1el.msh" );
   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                    "../../data/meshes/3D/pyramid_2el.msh" );
   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                    "../../data/meshes/3D/pyramid_2el.msh" );

   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                    "../../data/meshes/3D/pyramid_4el.msh" );
   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                    "../../data/meshes/3D/pyramid_4el.msh" );
   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI,
                                    "../../data/meshes/3D/regular_octahedron_8el.msh" );
   testN1E1AdditiveCommunication3D( communication::BufferedCommunicator::LocalCommunicationMode::DIRECT,
                                    "../../data/meshes/3D/regular_octahedron_8el.msh" );

   return EXIT_SUCCESS;
}
