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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/gridtransferoperators/Embeddings.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

static void testP2P1Transfer()
{
   const uint_t level = 4;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto p1Function = std::make_shared< P1Function< real_t > >( "p1Function", storage, level, level );
   auto p2Function = std::make_shared< P2Function< real_t > >( "p2Function", storage, level, level );

   VTKOutput vtkOutput( "../../output", "P2P1TransferTest", storage );
   vtkOutput.add( *p1Function );
   vtkOutput.add( *p2Function );

   // To test the transfer we
   // 1. set one vertex unknown in the middle of a P1 macro-face to testValue (rest 0.0)
   // 2. we prolongate to P2 and should get testValue at the corresponding vertex unknown and 0.5 * testValue at all (six, directly) neighboring edge unknowns
   // 3. we restrict again and should get a 2.5 * testValue at the corresponding vertex unknown and 0.25 * testValue at all neighboring vertex unknowns

   // Step 1:

   std::function< real_t( const Point3D& ) > zeros = []( const Point3D& ) { return 0; };
   p1Function->interpolate( zeros, level );
   p2Function->interpolate( zeros, level );

   const uint_t x         = 4;
   const uint_t y         = 4;
   const real_t testValue = 1.0;

   std::vector< PrimitiveID > faceIDs;
   storage->getFaceIDs( faceIDs );

   WALBERLA_CHECK_EQUAL( storage->getNumberOfLocalFaces(), 1 );
   WALBERLA_CHECK_EQUAL( faceIDs.size(), 1 );

   const auto p1FaceDataID = p1Function->getFaceDataID();
   auto       p1FaceData   = storage->getFace( faceIDs[0] )->getData( p1FaceDataID )->getPointer( level );

   const uint_t idx = vertexdof::macroface::indexFromVertex( level, x, y, stencilDirection::VERTEX_C );

   p1FaceData[idx] = testValue;

   vtkOutput.write( level, 1 );

   // Step 2:

   p2Function->prolongateP1ToP2( *p1Function, level );

   const auto p2VertexDoFFaceDataID = p2Function->getVertexDoFFunction().getFaceDataID();
   const auto p2EdgeDoFFaceDataID   = p2Function->getEdgeDoFFunction().getFaceDataID();

   auto p2VertexDoFFaceData = storage->getFace( faceIDs[0] )->getData( p2VertexDoFFaceDataID )->getPointer( level );
   auto p2EdgeDoFFaceData   = storage->getFace( faceIDs[0] )->getData( p2EdgeDoFFaceDataID )->getPointer( level );

   WALBERLA_CHECK_FLOAT_EQUAL( p2VertexDoFFaceData[idx], testValue )

   WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_W )],
                               0.5 * testValue )
   WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_E )],
                               0.5 * testValue )
   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_NW )], 0.5 * testValue )
   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_SE )], 0.5 * testValue )
   WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_N )],
                               0.5 * testValue )
   WALBERLA_CHECK_FLOAT_EQUAL( p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_S )],
                               0.5 * testValue )

   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_NW )], 0.0 )
   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_HO_SE )], 0.0 )
   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_SW )], 0.0 )
   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_DI_NE )], 0.0 )
   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_NW )], 0.0 )
   WALBERLA_CHECK_FLOAT_EQUAL(
       p2EdgeDoFFaceData[edgedof::macroface::indexFromVertex( level, x, y, stencilDirection::EDGE_VE_SE )], 0.0 )

   vtkOutput.write( level, 2 );

   // Step 3:

   p2Function->restrictP2ToP1( *p1Function, level );

   WALBERLA_CHECK_FLOAT_EQUAL( p1FaceData[idx], 2.5 * testValue );

   for ( const auto& neighbor : vertexdof::macroface::neighborsWithoutCenter )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( p1FaceData[vertexdof::macroface::indexFromVertex( level, x, y, neighbor )], 0.25 * testValue );
   }

   vtkOutput.write( level, 3 );
}

} // namespace hyteg

using namespace hyteg;

void logSectionHeader( const char* header )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separator( len + 2, '-' );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << "\n" << separator );
}

void run_P1ToP2EmbeddingTest( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   // auto linearPolynomial = []( const Point3D& x ) { return real_c( 3 ) * x[0] + real_c( 0.5 ) * x[1] - real_c( 1 ); };
   auto linearPolynomial = []( const Point3D& ) { return real_c( 1 ); };

   P1Function< real_t > p1Func( "p1", storage, level, level );
   P2Function< real_t > p2FuncDirect( "direct", storage, level, level );
   P2Function< real_t > p2FuncEmbedded( "via P1", storage, level, level );
   P2Function< real_t > p2FuncError( "error", storage, level, level );

   p1Func.interpolate( linearPolynomial, level, All );
   p2FuncDirect.interpolate( linearPolynomial, level, All );

   embedP1IntoP2( p1Func, p2FuncEmbedded, level, All );

   p2FuncError.assign( { real_c( 1 ), real_c( -1 ) }, { p2FuncDirect, p2FuncEmbedded }, level, All );

   bool outputVTK{ true };

   if ( outputVTK )
   {
      VTKOutput vtkOutput( "../../output", "P1ToP2EmbeddingOperatorTest", storage );
      vtkOutput.add( p1Func );
      vtkOutput.add( p2FuncDirect );
      vtkOutput.add( p2FuncEmbedded );
      vtkOutput.add( p2FuncError );
      vtkOutput.write( level );
   }

   real_t errorMeasure = p2FuncError.getMaxMagnitude( level, All );
   WALBERLA_LOG_INFO_ON_ROOT( " measure of error = " << errorMeasure );
   real_t tol{ std::is_same< double, real_t >() ? 1e-12 : 1e-6f };
   WALBERLA_CHECK_LESS( errorMeasure, tol );
}

void run2D_P1ToP2EmbeddingTest( uint_t level )
{
   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, -1.0 ), Point2D( 2.0, 3.0 ), MeshInfo::CRISS, 1, 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage );

   run_P1ToP2EmbeddingTest( storage, level );
}

void run3D_P1ToP2EmbeddingTest( uint_t level )
{
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/regular_octahedron_8el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage );

   run_P1ToP2EmbeddingTest( storage, level );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   logSectionHeader( "Testing transfers in 2D" );
   hyteg::testP2P1Transfer();

   logSectionHeader( "Testing embedding in 2D" );
   run2D_P1ToP2EmbeddingTest( 3 );

   // Commented out, see issue #205
   // logSectionHeader( "Testing embedding in 3D" );
   // run3D_P1ToP2EmbeddingTest( 2 );

   return EXIT_SUCCESS;
}
