/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Marcus Mohr.
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
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

// to be removed when moving to walberla namespace
using walberla::mpi::RecvBuffer;
using walberla::mpi::SendBuffer;

static void testPrimitiveID()
{
   PrimitiveID id0;
   PrimitiveID id1 = PrimitiveID::create( 42 );
   PrimitiveID id2( id1 );

   SendBuffer sendBuffer;

   PrimitiveID id3, id4, id5;

   sendBuffer << id0 << id1 << id2;

   RecvBuffer recvBuffer( sendBuffer );

   recvBuffer >> id3 >> id4 >> id5;

   WALBERLA_CHECK_EQUAL( id3, id0 );
   WALBERLA_CHECK_EQUAL( id4, id1 );
   WALBERLA_CHECK_EQUAL( id5, id2 );

   auto parent = id1;
   for ( uint_t k = 1; k <= PrimitiveID::maxRefinement(); ++k )
   {
      auto children = parent.createChildren();

      SendBuffer snd;
      for ( auto& child : children )
      {
         WALBERLA_CHECK_EQUAL( child.numAncestors(), k );
         WALBERLA_CHECK_EQUAL( child.getParent(), parent );
         snd << child;
      }

      RecvBuffer rcv( snd );
      for ( auto& child : children )
      {
         rcv >> id5;
         WALBERLA_CHECK_EQUAL( child, id5 );
      }

      parent = children[0];
   }
}

static void testPrimitiveTypes()
{
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
   PrimitiveStorage storage( setupStorage );

   for ( const auto& it : storage.getVertices() )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Primitive with ID = " << it.first << " is a vertex" );
      WALBERLA_CHECK_EQUAL( it.second->getType(), Primitive::PrimitiveTypeEnum::VERTEX );
   }
   for ( const auto& it : storage.getEdges() )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Primitive with ID = " << it.first << " is an edge" );
      WALBERLA_CHECK_EQUAL( it.second->getType(), Primitive::PrimitiveTypeEnum::EDGE );
   }
   for ( const auto& it : storage.getFaces() )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Primitive with ID = " << it.first << " is a face" );
      WALBERLA_CHECK_EQUAL( it.second->getType(), Primitive::PrimitiveTypeEnum::FACE );
   }
   for ( const auto& it : storage.getCells() )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Primitive with ID = " << it.first << " is a cell" );
      WALBERLA_CHECK_EQUAL( it.second->getType(), Primitive::PrimitiveTypeEnum::CELL );
   }
}

void testPrimitiveIDFormatting()
{
   // Generate a small PrimitiveStorage
   MeshInfo              meshInfo = MeshInfo::singleTriangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 0.0 ), Point2D( 1.0, 1.0 ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

   std::stringstream sstr;

   std::string control{
       "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000001034500000000000000000000000000000110" };

   for ( const auto& item : storage->getVertices() )
   {
      WALBERLA_LOG_INFO( "PrimitiveStorage holds VERTEX " << item.second->getID() );
      sstr << item.second->getID();
   }

   for ( const auto& item : storage->getEdges() )
   {
      WALBERLA_LOG_INFO( "PrimitiveStorage holds EDGE   " << PrimitiveIDFormat::COARSE_ID << item.second->getID() );
      sstr << item.second->getID();
   }

   PrimitiveIDFormatter::set( PrimitiveIDFormat::TRIM );
   for ( const auto& item : storage->getFaces() )
   {
      WALBERLA_LOG_INFO( "PrimitiveStorage holds FACE   " << item.second->getID() );
      sstr << item.second->getID();
   }

   WALBERLA_CHECK_EQUAL( sstr.str(), control );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "Running testPrimitiveID()" );
   hyteg::testPrimitiveID();

   WALBERLA_LOG_INFO_ON_ROOT( "Running testPrimitiveTypes()" );
   hyteg::testPrimitiveTypes();

   WALBERLA_LOG_INFO_ON_ROOT( "Running testPrimitiveIDFormatting()" );
   hyteg::testPrimitiveIDFormatting();

   return EXIT_SUCCESS;
}
