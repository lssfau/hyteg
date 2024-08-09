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
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

class TestData
{
 public:
   bool                a = false;
   int                 i = 100;
   std::vector< bool > aa;
};

class VertexTestData
{
 public:
   ~VertexTestData() { delete[] f; }

   bool                a = false;
   int                 i = 200;
   std::vector< bool > aa;
   float*              f;
};

class EdgeTestData
{
 public:
   bool                a = false;
   int                 i = 300;
   std::vector< bool > aa;
};

class CellTestData
{
 public:
   bool                a = false;
   int                 i = 400;
   std::vector< bool > aa;
};

class TestDataHandling : public OnlyInitializeDataHandling< TestData, Primitive >
{
 public:
   std::shared_ptr< TestData > initialize( const Primitive* const ) const
   {
      auto testData = std::make_shared< TestData >();
      testData->i   = 5555;
      return testData;
   }
};

class VertexTestDataHandling : public OnlyInitializeDataHandling< VertexTestData, Vertex >
{
 public:
   std::shared_ptr< VertexTestData > initialize( const Vertex* const ) const
   {
      auto testData = std::make_shared< VertexTestData >();
      testData->i   = 6666;
      testData->f   = new float[10000];
      return testData;
   }
};

class EdgeTestDataHandling : public OnlyInitializeDataHandling< EdgeTestData, Edge >
{
 public:
   std::shared_ptr< EdgeTestData > initialize( const Edge* const ) const
   {
      auto testData = std::make_shared< EdgeTestData >();
      testData->i   = 7777;
      return testData;
   }
};

class CellTestDataHandling : public OnlyInitializeDataHandling< CellTestData, Cell >
{
 public:
   std::shared_ptr< CellTestData > initialize( const Cell* const ) const
   {
      auto testData = std::make_shared< CellTestData >();
      testData->i   = 9999;
      return testData;
   }
};

static void testPrimitiveData()
{
   const std::string     meshFileName = prependHyTeGMeshDir( "3D/cube_24el.msh" );
   MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   WALBERLA_LOG_INFO( setupStorage );

   PrimitiveStorage storage( setupStorage );

   auto testDataHandling       = std::make_shared< TestDataHandling >();
   auto vertexTestDataHandling = std::make_shared< VertexTestDataHandling >();
   auto edgeTestDataHandling   = std::make_shared< EdgeTestDataHandling >();
   auto cellTestDataHandling   = std::make_shared< CellTestDataHandling >();

   // Adding data to all primitives
   PrimitiveDataID< TestData, Primitive > testDataID;
   storage.addPrimitiveData( testDataID, testDataHandling, "primitive data" );
   // Adding data only to vertices
   PrimitiveDataID< VertexTestData, Vertex > vertexTestDataID;
   storage.addVertexData( vertexTestDataID, vertexTestDataHandling, "vertex data" );
   // Adding data only to cells
   PrimitiveDataID< CellTestData, Cell > cellTestDataID;
   storage.addCellData( cellTestDataID, cellTestDataHandling, "cell data" );

   std::vector< PrimitiveID > primitiveIDs;
   storage.getPrimitiveIDs( primitiveIDs );
   for ( const auto& id : primitiveIDs )
   {
      WALBERLA_LOG_PROGRESS( "Checking content of primitive with ID: " << id );
      auto      primitive = storage.getPrimitive( id );
      TestData* testData  = primitive->getData( testDataID );
      WALBERLA_CHECK_EQUAL( testData->i, 5555 );
   }

   // Obtaining initialized vertex data from a vertex
   for ( const auto& it : storage.getVertices() )
   {
      WALBERLA_LOG_PROGRESS( "Checking content of vertex with ID: " << it.second->getID() );
      VertexTestData* vertexTestData = it.second->getData( vertexTestDataID );
      WALBERLA_CHECK_EQUAL( vertexTestData->i, 6666 );
   }

   // Obtaining initialized cell data from a cell
   for ( const auto& it : storage.getCells() )
   {
      WALBERLA_LOG_PROGRESS( "Checking content of cell with ID: " << it.second->getID() );
      CellTestData* cellTestData = it.second->getData( cellTestDataID );
      WALBERLA_CHECK_EQUAL( cellTestData->i, 9999 );
   }

   WALBERLA_UNUSED( testDataID );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testPrimitiveData();

   return EXIT_SUCCESS;
}
