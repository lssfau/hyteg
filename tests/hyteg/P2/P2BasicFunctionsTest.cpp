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
#include <functional>
#include <vector>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/PrimitiveID.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;

namespace hyteg {

static void testP2BasicFunctions()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P2Function< real_t > x( "x", storage, minLevel, maxLevel );
   P2Function< real_t > y( "y", storage, minLevel, maxLevel );
   P2Function< real_t > z( "y", storage, minLevel, maxLevel );

   std::vector< PrimitiveID > faces;
   storage->getFaceIDs( faces );
   Face* face = storage->getFace( faces[0] );

   real_t* faceEdgeDataX = face->getData( x.getEdgeDoFFunction().getFaceDataID() )->getPointer( maxLevel );
   real_t* faceEdgeDataY = face->getData( y.getEdgeDoFFunction().getFaceDataID() )->getPointer( maxLevel );
   real_t* faceEdgeDataZ = face->getData( z.getEdgeDoFFunction().getFaceDataID() )->getPointer( maxLevel );

   real_t* faceVertexDataX = face->getData( x.getVertexDoFFunction().getFaceDataID() )->getPointer( maxLevel );
   real_t* faceVertexDataY = face->getData( y.getVertexDoFFunction().getFaceDataID() )->getPointer( maxLevel );
   real_t* faceVertexDataZ = face->getData( z.getVertexDoFFunction().getFaceDataID() )->getPointer( maxLevel );

   // Interpolate

   std::function< real_t( const hyteg::Point3D& ) > expr  = []( const Point3D& ) -> real_t { return real_c( 2 ); };
   std::function< real_t( const hyteg::Point3D& ) > zeros = []( const Point3D& ) -> real_t { return real_c( 0 ); };
   std::function< real_t( const hyteg::Point3D& ) > func  = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + xx[0] ) * ( 2.0 + xx[1] ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > func2 = []( const Point3D& xx ) -> real_t {
      return real_c( ( 1.0 + ( xx[0] / 5.0 ) ) * ( 2.0 + xx[1] ) );
   };

   walberla::WcTimingPool timer;

   timer["Interpolate"].start();
   x.interpolate( expr, maxLevel, DoFType::All );
   y.interpolate( expr, maxLevel, DoFType::All );
   z.interpolate( func, maxLevel, DoFType::All );
   timer["Interpolate"].end();

   hyteg::communication::syncP2FunctionBetweenPrimitives( x, maxLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataX[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataX[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )], real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataX[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )], real_c( 2 ) );

      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )], real_c( 2 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )], real_c( 2 ) );
   }

   hyteg::communication::syncP2FunctionBetweenPrimitives( x, maxLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataX[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 2 ) );

      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 2 ) );
   }

   // Assign

   timer["Assign"].start();
   y.assign( {3.0, 2.0}, {x, y}, maxLevel, DoFType::All );
   z.assign( {1.0, -1.0}, {z, z}, maxLevel, DoFType::All );
   timer["Assign"].end();

   hyteg::communication::syncP2FunctionBetweenPrimitives( y, maxLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( z, maxLevel );

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 10 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 10 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 10 ) );
   }

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 10 ) );
   }

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 0 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )], real_c( 0 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )], real_c( 0 ) );
   }

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataZ[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 0 ) );
   }

   // Add

   timer["Add"].start();
   y.add( {4.0, 3.0}, {x, y}, maxLevel, DoFType::All );
   timer["Add"].end();

   hyteg::communication::syncP2FunctionBetweenPrimitives( y, maxLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( z, maxLevel );

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 48 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 48 ) );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataY[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )],
                                  real_c( 48 ) );
   }

   // Dot

   timer["Dot"].start();
   const real_t scalarProduct = y.dotGlobal( x, maxLevel, DoFType::All );
   timer["Dot"].end();

   uint_t totalPointsOnFace = levelinfo::num_microvertices_per_face( maxLevel ) + levelinfo::num_microedges_per_face( maxLevel );
   WALBERLA_CHECK_FLOAT_EQUAL( scalarProduct, real_c( totalPointsOnFace * 48 * 2 ) );

   WALBERLA_LOG_INFO_ON_ROOT( timer );

   ///Check assign against add

   x.interpolate( func, maxLevel, DoFType::All );
   y.interpolate( func2, maxLevel, DoFType::All );
   z.interpolate( zeros, maxLevel, DoFType::All );

   z.assign( {1.0, -1.0}, {x, y}, maxLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( z, maxLevel );
   x.add( {-1.0}, {y}, maxLevel );
   hyteg::communication::syncP2FunctionBetweenPrimitives( x, maxLevel );

   for( const auto& it : edgedof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )],
                                  faceEdgeDataX[edgedof::macroface::horizontalIndex( maxLevel, it.x(), it.y() )] );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )],
                                  faceEdgeDataX[edgedof::macroface::diagonalIndex( maxLevel, it.x(), it.y() )] );
      WALBERLA_CHECK_FLOAT_EQUAL( faceEdgeDataZ[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )],
                                  faceEdgeDataX[edgedof::macroface::verticalIndex( maxLevel, it.x(), it.y() )] );
   }

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataZ[vertexdof::macroface::index( maxLevel, it.x(), it.y() )],
                                  faceVertexDataX[vertexdof::macroface::index( maxLevel, it.x(), it.y() )] );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2BasicFunctions();

   return EXIT_SUCCESS;
}
