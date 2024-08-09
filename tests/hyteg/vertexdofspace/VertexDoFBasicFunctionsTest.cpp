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
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_c;

static void testVertexDoFBasicFunctions()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   MeshInfo mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto x = vertexdof::VertexDoFFunction< real_t >( "x", storage, minLevel, maxLevel );
   auto y = vertexdof::VertexDoFFunction< real_t >( "y", storage, minLevel, maxLevel );

   std::vector< PrimitiveID > faces;
   storage->getFaceIDs( faces );
   Face* face = storage->getFace( faces[0] );

   real_t* faceVertexDataX = face->getData( x.getFaceDataID() )->getPointer( maxLevel );
   real_t* faceVertexDataY = face->getData( y.getFaceDataID() )->getPointer( maxLevel );

   // Interpolate

   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& ) -> real_t { return real_c( 2 ); };

   walberla::WcTimingPool timer;

   timer["Interpolate"].start();
   x.interpolate( expr, maxLevel, DoFType::All );
   y.interpolate( expr, maxLevel, DoFType::All );
   timer["Interpolate"].end();

   hyteg::communication::syncFunctionBetweenPrimitives( x, maxLevel );
   hyteg::communication::syncFunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataX[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 2 ) );

      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 2 ) );
   }

   // Assign

   timer["Assign"].start();
   //y->assign( { 3.0, 2.0 }, { x.get(), y.get() }, maxLevel, DoFType::All );
   y.assign( {{3.0, 2.0}}, {{x, y}}, maxLevel, DoFType::All );
   timer["Assign"].end();

   hyteg::communication::syncFunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 10 ) );
   }

   // Add

   timer["Add"].start();
   y.add( {{4.0, 3.0}}, {{x, y}}, maxLevel, DoFType::All );
   timer["Add"].end();

   hyteg::communication::syncFunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 48 ) );
   }

   // Dot

   timer["Dot"].start();
   const real_t scalarProduct = y.dotGlobal( x, maxLevel, DoFType::All );
   timer["Dot"].end();

   WALBERLA_CHECK_FLOAT_EQUAL( scalarProduct, real_c( levelinfo::num_microvertices_per_face( maxLevel ) * 48 * 2 ) );

   timer["MultElementWise"].start();
   y.multElementwise( {{x, y}}, maxLevel, DoFType::All );
   timer["MultElementWise"].end();
   hyteg::communication::syncFunctionBetweenPrimitives( y, maxLevel );

   for( const auto& it : vertexdof::macroface::Iterator( maxLevel ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( faceVertexDataY[vertexdof::macroface::index( maxLevel, it.x(), it.y() )], real_c( 96 ) );
   }

   WALBERLA_LOG_INFO_ON_ROOT( timer );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testVertexDoFBasicFunctions();

   return EXIT_SUCCESS;
}
