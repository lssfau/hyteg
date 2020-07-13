/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#pragma once

#include <vector>

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

namespace hyteg {

using namespace walberla::mpistubs;
using walberla::int_c;
using walberla::uint_c;
using walberla::uint_t;
using walberla::mpi::MPIRank;

/// \brief BufferSystem alternative if not the sender ranks but the number of expected messages is known.
template < typename RecvBuffer_T = walberla::mpi::RecvBuffer, typename SendBuffer_T = walberla::mpi::SendBuffer >
class GenericPackageBufferSystem
{
 public:
   /// \brief Constructs a new PackageBufferSystem
   ///
   /// \param comm The MPI communicator to work on.
   /// \param numExpectedPackages The (exact) number of packages that this rank is expected to receive.
   ///                            The sender rank is irrelevant in this case.
   /// \param tag A tag for the MPI communication.
   ///
   GenericPackageBufferSystem( const MPI_Comm& comm, uint_t numExpectedPackages, int tag = 0 )
   : comm_( comm )
   , numExpectedPackages_( numExpectedPackages )
   , tag_( tag )
   , numPackagesReceived_( 0 )
   , allPackagesSent_( false )
   {}

   class PackageRecvInfo
   {
    public:
      friend class GenericPackageBufferSystem;

      RecvBuffer_T& buffer() { return buffer_; }
      uint_t        size() const { return size_; }
      MPIRank       senderRank() const { return senderRank_; }

    private:
      RecvBuffer_T buffer_;
      uint_t       size_;
      MPIRank      senderRank_;
   };

   /// \brief Returns the buffer of a new package that shall be sent.
   ///
   /// Each call to this function creates a new buffer that represents
   /// one package that is then received by the specified target rank.
   ///
   /// \param targetRank The receiver rank of this package.
   /// \return Reference to a send buffer associated with a package.
   ///         Be careful to store this buffer in a reference if not used directly.
   ///
   SendBuffer_T& getPackageSendBuffer( MPIRank targetRank )
   {
      WALBERLA_CHECK( !allPackagesSent_,
                      "Reuse of PackageBufferSystem not supported. Create a new one or implement reset() functionality." )

      PackageSendInfo newSendInfo( targetRank );
      sendInfos_.push_back( newSendInfo );
      return sendInfos_.back().buffer;
   }

   SendBuffer_T& getPackageSendBuffer( uint_t targetRank ) { return getPackageSendBuffer( MPIRank( targetRank ) ); }

   /// \brief Sends all packages.
   ///
   /// This BufferSystem cannot be reused after this method has been called.
   /// Create a new instance if you want to perform another communication step.
   ///
   void sendAll()
   {
      allPackagesSent_ = true;

      WALBERLA_NON_MPI_SECTION()
      {
         return;
      }

      MPI_Request request;

      for ( const auto& it : sendInfos_ )
      {
         MPI_Isend( it.buffer.ptr(),           // pointer to size buffer
                    int_c( it.buffer.size() ), // send one size
                    MPI_BYTE,                  // type
                    it.targetRank,             // receiver rank
                    tag_,                      // message tag
                    comm_,                     // communicator
                    &request                   // request needed for wait
         );
      }
   }

   /// \brief Returns false as long as not all packages have been obtained via getNextPackage().
   bool allPackagesReceived() const { return numPackagesReceived_ >= numExpectedPackages_; }

   /// \brief Returns the next received package.
   ///
   /// \return A wrapper class that contains all relevant information about the package as well
   ///         as an associated receive buffer.
   PackageRecvInfo getNextPackage()
   {
      WALBERLA_CHECK( !allPackagesReceived() );

      WALBERLA_NON_MPI_SECTION()
      {
         PackageSendInfo& sendInfo   = sendInfos_[numPackagesReceived_];
         const auto       senderRank = walberla::mpi::MPIManager::instance()->rank();
         WALBERLA_CHECK_EQUAL( sendInfo.targetRank, senderRank );
         PackageRecvInfo recvInfo;

         recvInfo.buffer()    = RecvBuffer_T( sendInfo.buffer );
         recvInfo.senderRank_ = senderRank;
         recvInfo.size_       = sendInfo.buffer.size();
         numPackagesReceived_++;
         return recvInfo;
      }

      PackageRecvInfo recvInfo;

      while ( true )
      {
         int        probeFlag;
         MPI_Status probeStatus;
         MPI_Iprobe( MPI_ANY_SOURCE, tag_, comm_, &probeFlag, &probeStatus );
         if ( probeFlag )
         {
            MPIRank sender = probeStatus.MPI_SOURCE;
            int     count  = 0;
            MPI_Get_count( &probeStatus, MPI_BYTE, &count );
            recvInfo.senderRank_ = sender;
            recvInfo.size_       = uint_c( count );
            recvInfo.buffer_.resize( recvInfo.size() );

            MPI_Status recvStatus;
            MPI_Recv( recvInfo.buffer().ptr(), // where to store received size
                      count,                   // size of expected message
                      MPI_BYTE,                // type
                      sender,                  // rank of sender process
                      tag_,                    // message tag
                      comm_,                   // communicator
                      &recvStatus              // request, needed for wait
            );

            numPackagesReceived_++;
            return recvInfo;
         }
      }
   }

 private:
   struct PackageSendInfo
   {
      PackageSendInfo( MPIRank _targetRank )
      : targetRank( _targetRank )
      {}
      MPIRank      targetRank;
      SendBuffer_T buffer;
   };

   MPI_Comm comm_;
   uint_t   numExpectedPackages_;
   int      tag_;

   std::vector< PackageSendInfo > sendInfos_;
   uint_t                         numPackagesReceived_;
   bool                           allPackagesSent_;
};

typedef GenericPackageBufferSystem<> PackageBufferSystem;

} // namespace hyteg