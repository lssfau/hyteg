/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/communication/BufferedCommunication.hpp"
#include "hyteg/communication/PackInfo.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::uint_t;

namespace hyteg {

/**
 * \page BH.02_Communication Tutorial BH.02 - Communicating data between primitives
 *
 * \dontinclude tutorials/basics-of-hyteg/BH.02_Communication.cpp
 *
 * \brief In this tutorial we will communicate data between primitives
 *
 * \section BH02-Communication-intro Introduction
 *
 * During typical simulations, the simulation domain is partitioned and distributed
 * among processes. In hhg, the primitives define the topology of the domain.
 *
 * Each primitive is assigned to one process, while one process may carry more
 * than one primitive.
 *
 * Since most simulations require data to be exchanged between primitives
 * (think of the ghost layer / halo concept) we need to be able to support
 * communcation between primitives that reside on different or the same process.
 *
 * The general approach can be divided into three steps:
 * - packing the data into send buffers on the sender side
 * - sending the buffer to a receiving primitive
 * - unpacking the data on the receiver side.
 *
 * We will in this tutorial show how to implement these steps.
 *
 * \section BH02-Communication-packinfo (Un-)Packing data
 *
 * For packing and unpacking data into and from MPI buffers, the framework provides the
 * interface hyteg::communication::PackInfo.
 *
 * All virtual methods must be implemented for a specific data item that shall be
 * transferred among different primitive types.
 *
 * Note that by design, it is not allowed to transfer data between two primitives
 * of the same type (e.g. vertex to edge is allowed while vertex to vertex is not).
 *
 * To keep things simple let's on every edge collect the IDs of the neighboring vertices via communication
 * (this is of course not necessary in real applications but gives a simple example).
 *
 * Therefore we introduce a struct TestData (plus the respective data handling) that carries the ID of the
 * carrier and a vector of neighboring IDs:
 *
 * \snippet{trimleft} this TestData
 *
 * Now let us implement the three virtual functions of the abstract class hyteg::communication::PackInfo that are needed to
 * send data from vertices to edges:
 *
 * \snippet{trimleft} this PackInfo
 *
 * Note that we implement the method communicateLocalVertexToEdge. It can be used automatically (we will see that later)
 * if both primitives reside on the same process to optimize the communication step.
 *
 * \section BH02-Communication Communication 
 *
 * To communicate the data, we need a hyteg::PrimitiveStorage and primitives that carry our struct:
 *
 * \snippet{trimleft} this Setup
 *
 * Then we create an instance of our PackInfo implementation and an instance of hyteg::communication::BufferedCommunicator
 * which will carry out the communication.
 *
 * \snippet{trimleft} this Communicator
 *
 * To tell the communicator that it shall communicate our test data during every communication, we need
 * to add our PackInfo. Only data with registered PackInfo will be communicated.
 *
 * \snippet{trimleft} this AddPackInfo
 *
 * As previously mentioned, local communication can be performed directly, without buffering (which is set by default).
 * However, the local communication mode can be changed via:
 *
 * \snippet{trimleft} this LocalMode
 *
 * See hyteg::communication::BufferedCommunicator for more details on the local communication modes.
 *
 * Now we perform the buffered, non-blocking communication:
 *
 * \snippet{trimleft} this Communication
 *
 * Since the communication is buffered and non-blocking, we perform it in two steps:
 * - startCommunication packs the data into buffers and starts non-blocking MPI sends
 * - endCommunication waits for the MPI receives and unpacks the data
 *
 * Therefore we can right after startCommunication() safely operate on the data of the sender side,
 * in our case on the vertex data. This allows for overlapping communication and computation.
 *
 * After the communication is done, we expect to have two IDs unpacked on all edges since every
 * edge has two neighboring vertices. Let's check that:
 *
 * \snippet{trimleft} this Check
 *
 * \section BH02-Communication-code Complete Program
 *
 * \include tutorials/basics-of-hyteg/BH.02_Communication.cpp
 *
 */

/// [TestData]
struct TestData
{
   PrimitiveID                ownID;
   std::vector< PrimitiveID > neighborIDs;
};

struct TestDataHandling : OnlyInitializeDataHandling< TestData, Primitive >
{
   virtual std::shared_ptr< TestData > initialize( const Primitive* const primitive ) const
   {
      auto data   = std::make_shared< TestData >();
      data->ownID = primitive->getID();
      return data;
   }
};
/// [TestData]

/// [PackInfo]
class TestPackInfo : public communication::PackInfo
{
 public:
   TestPackInfo( PrimitiveDataID< TestData, Primitive >& dataID )
   : dataID_( dataID )
   {}

   virtual void
       packVertexForEdge( const Vertex* sender, const PrimitiveID& /* receiver */, walberla::mpi::SendBuffer& buffer ) const
   {
      WALBERLA_LOG_INFO( "Packing data on vertex..." );

      TestData* data = sender->getData( dataID_ );

      // The operator<< overload allows for easy buffer packing
      // of standard data types (thanks to the waLBerla framework).
      buffer << data->ownID;
   }

   virtual void unpackEdgeFromVertex( Edge* receiver, const PrimitiveID& /* sender */, walberla::mpi::RecvBuffer& buffer ) const
   {
      WALBERLA_LOG_INFO( "Unpacking data on edge..." );

      TestData*   data = receiver->getData( dataID_ );
      PrimitiveID vertexData;
      buffer >> vertexData;

      // Adding the received ID to the neighbors..
      data->neighborIDs.push_back( vertexData );
   }

   virtual void communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
   {
      WALBERLA_LOG_INFO( "Communicating data unbuffered from vertex to edge..." );

      TestData* vertexData = sender->getData( dataID_ );
      TestData* edgeData   = receiver->getData( dataID_ );
      edgeData->neighborIDs.push_back( vertexData->ownID );
   }

   // Left other methods empty for this tutorial.

 private:
   PrimitiveDataID< TestData, Primitive > dataID_;

   /// [PackInfo]

 public:
   virtual void packEdgeForVertex( const Edge*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const {}

   virtual void unpackVertexFromEdge( Vertex*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const {}

   virtual void communicateLocalEdgeToVertex( const Edge*, Vertex* ) const {}

   virtual void packEdgeForFace( const Edge*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const {}

   virtual void unpackFaceFromEdge( Face*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const {}

   virtual void communicateLocalEdgeToFace( const Edge*, Face* ) const {}

   virtual void packFaceForEdge( const Face*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const {}

   virtual void unpackEdgeFromFace( Edge*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const {}

   virtual void communicateLocalFaceToEdge( const Face*, Edge* ) const {}
};

void CommunicationTutorial()
{
   uint_t numProcesses = walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "tri_2el.msh" ) );
   SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );

   loadbalancing::roundRobin( setupStorage );

   /// [Setup]
   // As seen in previous tutorials...
   auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   PrimitiveDataID< TestData, Primitive > dataID;
   auto                                   testDataHandling = std::make_shared< TestDataHandling >();
   storage->addPrimitiveData( dataID, testDataHandling, "test data" );
   /// [Setup]

   for ( const auto& it : storage->getEdges() )
   {
      auto edge = it.second;
      WALBERLA_CHECK_EQUAL( edge->getData( dataID )->neighborIDs.size(), 0 );
   }

   /// [Communicator]
   std::shared_ptr< TestPackInfo > packInfo = std::make_shared< TestPackInfo >( dataID );

   communication::BufferedCommunicator communicator( storage );
   /// [Communicator]

   /// [AddPackInfo]
   communicator.addPackInfo( packInfo );
   /// [AddPackInfo]

   /// [LocalMode]
   // communicator.setLocalCommunicationMode( communication::BufferedCommunicator::BUFFERED_MPI );
   /// [LocalMode]

   /// [Communication]
   // Communicate data from all vertices to the neighboring edges:

   // Packing data from vertices into buffers and starting non-blocking MPI calls
   communicator.startCommunication< Vertex, Edge >();

   // Waiting for MPI sends and unpacking data to all receiving edges
   communicator.endCommunication< Vertex, Edge >();
   /// [Communication]

   /// [Check]
   for ( const auto& it : storage->getEdges() )
   {
      auto edge = it.second;
      WALBERLA_CHECK_EQUAL( edge->getData( dataID )->neighborIDs.size(), 2 );
   }
   /// [Check]
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   hyteg::CommunicationTutorial();
   return 0;
}
