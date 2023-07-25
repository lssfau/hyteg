/*
 * Copyright (c) 2021-2023 Marcus Mohr
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
#include "hyteg/functions/FunctionMultiStore.hpp"

#include <vector>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::uint_t;

using namespace hyteg;

// --------------------------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   uint_t theLevel = 3;

   // Generate mesh around origin
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   meshInfo          = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CROSS, 1, 1 );

   // Generate primitives
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // =============
   //  P1Function
   // =============
   hyteg::P1Function< double > dFunc1( "#1", storage, theLevel, theLevel );
   hyteg::P1Function< double > dFunc2( "#2", storage, theLevel, theLevel );
   // hyteg::P1Function< float > fFunc1( "#3", storage, theLevel, theLevel );
   hyteg::P1Function< int32_t > iFunc1( "#4", storage, theLevel, theLevel );
   hyteg::P1Function< int64_t > lFunc1( "#5", storage, theLevel, theLevel );

   FunctionMultiStore< P1Function > ms;

   ms.push_back( dFunc1 );
   ms.push_back( dFunc2 );
   ms.push_back( iFunc1 );
   ms.push_back( lFunc1 );

   WALBERLA_LOG_INFO_ON_ROOT( "FunctionMultiStore holds " << ms.size() << " P1Functions" );
   WALBERLA_CHECK_EQUAL( ms.size(), 4 );

   // Remove a single of the P1Functions
   WALBERLA_LOG_INFO_ON_ROOT( "Removing '" << dFunc2.getFunctionName() << "' from FunctionMultiStore" );
   ms.remove( dFunc2 );
   WALBERLA_LOG_INFO_ON_ROOT( "FunctionMultiStore holds " << ms.size() << " P1Functions" );
   WALBERLA_CHECK_EQUAL( ms.size(), 3 );

   hyteg::P1Function< double > aux = dFunc1;

   // ==================
   //  P2VectorFunction
   // ==================
   hyteg::P2VectorFunction< double >  dVecFunc1( "#1", storage, theLevel, theLevel );
   hyteg::P2VectorFunction< double >  dVecFunc2( "#2", storage, theLevel, theLevel );
   hyteg::P2VectorFunction< double >  dVecFunc3( "#7", storage, theLevel, theLevel );
   hyteg::P2VectorFunction< float >   fVecFunc1( "#3", storage, theLevel, theLevel );
   hyteg::P2VectorFunction< int32_t > iVecFunc1( "#4", storage, theLevel, theLevel );
   hyteg::P2VectorFunction< int64_t > lVecFunc1( "#5", storage, theLevel, theLevel );
   hyteg::P2VectorFunction< int64_t > lVecFunc2( "#6", storage, theLevel, theLevel );

   FunctionMultiStore< P2VectorFunction > ms2;

   ms2.push_back( iVecFunc1 );
   ms2.push_back( lVecFunc1 );
   ms2.push_back( dVecFunc1 );
   ms2.push_back( dVecFunc2 );
   ms2.push_back( lVecFunc2 );
   ms2.push_back( dVecFunc3 );
   ms2.push_back( fVecFunc1 );

   WALBERLA_LOG_INFO_ON_ROOT( "FunctionMultiStore holds " << ms2.size() << " P2VectorFunctions" );

   // try obtaining a specific vector from the store
   auto dFuncs = ms2.getFunctions< double >();
   WALBERLA_LOG_INFO_ON_ROOT( " - of these " << dFuncs.size() << " have value type 'double'" );
   WALBERLA_CHECK_EQUAL( dFuncs.size(), 3 );

   auto iFuncs = ms2.getFunctions< int32_t >();
   WALBERLA_LOG_INFO_ON_ROOT( " - of these " << iFuncs.size() << " have value type 'int32_t'" );
   WALBERLA_CHECK_EQUAL( iFuncs.size(), 1 );

   auto lFuncs = ms2.getFunctions< int64_t >();
   WALBERLA_LOG_INFO_ON_ROOT( " - of these " << lFuncs.size() << " have value type 'int64_t'" );
   WALBERLA_CHECK_EQUAL( lFuncs.size(), 2 );

   auto fFuncs = ms2.getFunctions< float >();
   WALBERLA_LOG_INFO_ON_ROOT( " - of these " << fFuncs.size() << " have value type 'float'" );
   WALBERLA_CHECK_EQUAL( fFuncs.size(), 1 );

   WALBERLA_LOG_INFO_ON_ROOT( "Their names are:" );
   std::vector< std::string > names = ms2.getFunctionNames();
   std::sort( names.begin(), names.end() );
   uint_t count{0};
   for( const auto& name: names ) {
     WALBERLA_LOG_INFO_ON_ROOT( " - " << name );
     std::stringstream stream;
     stream << "#" << ++count;
     WALBERLA_CHECK_EQUAL( name, stream.str() );
   }

   // Empty the store of P2VectorFunctions
   ms2.remove( iVecFunc1 );
   WALBERLA_CHECK_EQUAL( ms2.size(), 6 );
   ms2.remove( lVecFunc1 );
   WALBERLA_CHECK_EQUAL( ms2.size(), 5 );
   ms2.remove( dVecFunc1 );
   WALBERLA_CHECK_EQUAL( ms2.size(), 4 );
   ms2.remove( dVecFunc2 );
   WALBERLA_CHECK_EQUAL( ms2.size(), 3 );
   ms2.remove( lVecFunc2 );
   WALBERLA_CHECK_EQUAL( ms2.size(), 2 );
   ms2.remove( dVecFunc3 );
   WALBERLA_CHECK_EQUAL( ms2.size(), 1 );
   ms2.remove( fVecFunc1 );
   WALBERLA_CHECK_EQUAL( ms2.size(), 0 );

   return EXIT_SUCCESS;
}
