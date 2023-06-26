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
#include "hyteg/p2functionspace/P2Function.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

static void testP2Function()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   MeshInfo mesh  = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   MeshInfo mesh2 = MeshInfo::fromGmshFile( "../../data/meshes/annulus_coarse.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   SetupPrimitiveStorage setupStorage2( mesh2, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage  = std::make_shared< PrimitiveStorage >( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage2 = std::make_shared< PrimitiveStorage >( setupStorage2 );

   auto x = P2Function< real_t >( "x", storage, minLevel, maxLevel );
   auto y = P2Function< real_t >( "y", storage, minLevel, maxLevel );

   std::vector< PrimitiveID > faces;
   storage->getFaceIDs( faces );
   Face* face = storage->getFace( faces[0] );

   real_t* vertexDoFFaceDataX = face->getData( x.getVertexDoFFunction().getFaceDataID() )->getPointer( maxLevel );
   real_t* vertexDoFFaceDataY = face->getData( y.getVertexDoFFunction().getFaceDataID() )->getPointer( maxLevel );
   real_t* edgeDoFFaceDataX   = face->getData( x.getEdgeDoFFunction().getFaceDataID() )->getPointer( maxLevel );
   real_t* edgeDoFFaceDataY   = face->getData( y.getEdgeDoFFunction().getFaceDataID() )->getPointer( maxLevel );

   // Interpolate

   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& ) -> real_t { return real_c( 2 ); };

   walberla::WcTimingPool timer;

   timer["Interpolate"].start();
   x.interpolate( expr, maxLevel, DoFType::All );
   timer["Interpolate"].end();

   hyteg::communication::syncP2FunctionBetweenPrimitives( x, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL(
          vertexDoFFaceDataX[vertexdof::macroface::indexFromVertex( maxLevel, it.x(), it.y(), stencilDirection::VERTEX_C )],
          real_c( 2 ) );
   }

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataX[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataX[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataX[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 2 ) );
   }

   // Assign

   timer["Assign"].start();
   y.assign( {3.0}, {x}, maxLevel, DoFType::All );
   timer["Assign"].end();

   hyteg::communication::syncP2FunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL(
          vertexDoFFaceDataY[vertexdof::macroface::indexFromVertex( maxLevel, it.x(), it.y(), stencilDirection::VERTEX_C )],
          real_c( 6 ) );
   }

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataY[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 6 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataY[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 6 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataY[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 6 ) );
   }

   // Add

   timer["Add"].start();
   y.add( {{4.0, 3.0}}, {{x, x}}, maxLevel, DoFType::All );
   timer["Add"].end();
   hyteg::communication::syncP2FunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL(
          vertexDoFFaceDataY[vertexdof::macroface::indexFromVertex( maxLevel, it.x(), it.y(), stencilDirection::VERTEX_C )],
          real_c( 20 ) );
   }

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataY[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 20 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataY[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 20 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFFaceDataY[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 20 ) );
   }

   // Dot

   timer["Dot"].start();
   const real_t scalarProduct = y.dotGlobal( x, maxLevel, DoFType::All );
   timer["Dot"].end();

   WALBERLA_CHECK_FLOAT_EQUAL(
       scalarProduct,
       real_c( ( levelinfo::num_microvertices_per_face( maxLevel ) + levelinfo::num_microedges_per_face( maxLevel ) ) * 20 *
               2 ) );

   WALBERLA_LOG_INFO_ON_ROOT( timer );

   // Output interpolate VTK

   auto p2 = std::make_shared< P2Function< real_t > >( "p2", storage2, minLevel, maxLevel );
   std::function< real_t( const hyteg::Point3D& ) > linearX = []( const Point3D& xx ) -> real_t { return xx[0] + xx[1]; };
   p2->interpolate( linearX, maxLevel, DoFType::All );

//   VTKOutput vtkOutput("../../output", "p2_interpolate_test", storage);
//   vtkOutput.add( *p2 );
//   vtkOutput.write( maxLevel );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2Function();

   return EXIT_SUCCESS;
}
