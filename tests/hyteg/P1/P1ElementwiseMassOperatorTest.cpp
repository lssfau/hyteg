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
// test that the product one^T*M*one with mass matrix M and vector of ones gives area of domain
#include "core/Environment.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/p1functionspace/P1ElementwiseOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;
using namespace hyteg;

void checkArea( std::shared_ptr< PrimitiveStorage > storage, real_t area )
{
   const size_t minLevel = 2;
   const size_t maxLevel = 4;

   hyteg::P1Function< real_t > microCoordX( "microCoordX", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > microCoordY( "microCoordY", storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > compX = []( const hyteg::Point3D& pp ) { return pp[0]; };
   std::function< real_t( const hyteg::Point3D& ) > compY = []( const hyteg::Point3D& pp ) { return pp[1]; };

   for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      microCoordX.interpolate( compX, lvl );
      microCoordY.interpolate( compY, lvl );

      communication::syncFunctionBetweenPrimitives( microCoordX, lvl );
      communication::syncFunctionBetweenPrimitives( microCoordY, lvl );
   }

   P1ElementwiseMassOperator massOp( storage, {&microCoordX, &microCoordY}, minLevel, maxLevel );

   P1Function< real_t >                      aux( "aux", storage, minLevel, maxLevel );
   P1Function< real_t >                      vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );
   std::function< real_t( const Point3D& ) > ones = []( const Point3D& ) { return 1.0; };

   for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      vecOfOnes.interpolate( ones, lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dotGlobal( aux, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "level " << lvl << ": measure = " << std::scientific << measure );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, area );
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Test with rectangle
   WALBERLA_LOG_INFO_ON_ROOT( "Testing with RECTANGLE" );
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0.0, 0.0} ), Point2D( {2.0, 1.0} ), MeshInfo::CRISS, 1, 1 );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   checkArea( storage, 2.0 );

   // Test with backward facing step
   WALBERLA_LOG_INFO_ON_ROOT( "Testing with BFS" );
   meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   SetupPrimitiveStorage setupStorageBFS( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storageBFS = std::make_shared< PrimitiveStorage >( setupStorageBFS );
   checkArea( storageBFS, 1.75 );

   return EXIT_SUCCESS;
}
