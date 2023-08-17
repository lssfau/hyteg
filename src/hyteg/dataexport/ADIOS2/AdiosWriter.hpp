/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include <adios2.h>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP1.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP2.hpp"
#include "hyteg/dataexport/FEFunctionRegistry.hpp"
#include "hyteg/dataexport/FEFunctionWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKMeshWriter.hpp"

namespace hyteg {

// using walberla::real_c;
// using walberla::uint64_t;
// using walberla::uint_c;
using walberla::real_t;
using walberla::uint_t;

class AdiosWriterForP1;
class AdiosWriterForP2;

/// Class to export FEFunction data using the ADIOS2 library
///
/// This class allows to write FEFunction data to the filesystem using the ADIOS2 library.
/// Output is using the BP format (currently version BP4). The class currently only supports
/// the following types of functions:
/// - P1Function and P1VectorFunction
/// - P2Function and P2VectorFunction
/// - P2P1TaylorHoodFunction
///
/// which must be using real_t as their value type.
class AdiosWriter : public FEFunctionWriter< AdiosWriter >
{
 public:
#ifdef WALBERLA_BUILD_WITH_MPI
   /// \param filePath          Path to directory where the BP files are stored
   /// \param fileBaseName      Basename for output BP file (which is acutally a folder)
   /// \param storage           PrimitiveStorage associated with functions to export
   /// \param comm              MPI Communicator, defaults to the HyTeG standard communicator
   AdiosWriter( std::string                                filePath,
                std::string                                fileBaseName,
                const std::shared_ptr< PrimitiveStorage >& storage,
                MPI_Comm                                   comm = walberla::MPIManager::instance()->comm() )
   : AdiosWriter( filePath, fileBaseName, "", storage, comm )
   {}

   /// \param filePath          Path to directory where the BP files are stored
   /// \param fileBaseName      Basename for output BP file (which is acutally a folder)
   /// \param configFile        Name of a file in XML or YAML format with runtime configuration parameters for ADIOS2
   /// \param storage           PrimitiveStorage associated with functions to export
   /// \param comm              MPI Communicator, defaults to the HyTeG standard communicator
   AdiosWriter( std::string                                filePath,
                std::string                                fileBaseName,
                std::string                                configFile,
                const std::shared_ptr< PrimitiveStorage >& storage,
                MPI_Comm                                   comm = walberla::MPIManager::instance()->comm() )
   : storage_( storage )
   , filePath_( filePath )
   , fileBaseName_( fileBaseName )
   {
      // setup central ADIOS2 interface object
      if ( configFile.empty() )
      {
         adios_ = adios2::ADIOS( comm );
      }
      else
      {
         adios_ = adios2::ADIOS( configFile, comm );
      }
   }

#else

   /// \param filePath          Path to directory where the BP files are stored
   /// \param fileBaseName      Basename for output BP file (which is acutally a folder)
   /// \param storage           PrimitiveStorage associated with functions to export
   AdiosWriter( std::string                                filePath,
                std::string                                fileBaseName,
                const std::shared_ptr< PrimitiveStorage >& storage )
   : AdiosWriter( filePath, fileBaseName, "", storage )
   {}

   /// \param filePath          Path to directory where the BP files are stored
   /// \param fileBaseName      Basename for output BP file (which is acutally a folder)
   /// \param configFile        Name of a file in XML or YAML format with runtime configuration parameters for ADIOS2
   /// \param storage           PrimitiveStorage associated with functions to export
   AdiosWriter( std::string                                filePath,
                std::string                                fileBaseName,
                std::string                                configFile,
                const std::shared_ptr< PrimitiveStorage >& storage )
   : storage_( storage )
   , filePath_( filePath )
   , fileBaseName_( fileBaseName )
   {
      // setup central ADIOS2 interface object
      if ( configFile.empty() )
      {
         adios_ = adios2::ADIOS();
      }
      else
      {
         adios_ = adios2::ADIOS( configFile );
      }
   }

#endif

   /// Add an FE Function to became part of the next dataexport phase
   template < template < typename > class func_t, typename value_t >
   inline void add( const func_t< value_t >& function )
   {
      // Still some places where we would need to add extra code to also
      // be able to write uint32/64_t
      if constexpr ( !std::is_same_v< real_t, value_t > )
      {
         WALBERLA_ABORT( "AdiosWriter currently only supports exporting FEFunction using real_t as value type!" );
      }

      if constexpr ( !std::is_same_v< func_t< value_t >, P1Function< value_t > > &&
                     !std::is_same_v< func_t< value_t >, P2Function< value_t > > &&
                     !std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > &&
                     !std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > &&
                     !std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
      {
         WALBERLA_ABORT( "AdiosWriter only supports P1 and P2 type functions!" );
      }

      if ( p1Writers_.size() > 0 || p2Writers_.size() > 0 )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter class does not support adding functions after the first write.\n"
                                       << "--> Ignoring function '" << function.getFunctionName() << "'!" );
      }
      else
      {
         feFunctionRegistry_.add( function );
      }
   }

   /// Set parameter specified by string key to value specified by string value
   ///
   /// Currently not supported for AdiosWriter
   void setParameter( const std::string& key, const std::string& value )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter::setParameter() does not perform any action!" );
   }

   void write( const uint_t level, const uint_t timestep = 0 );

   /// Class that wraps an Adios span such that we can insert data with operator<<
   ///
   /// The blockStride template argument is designed to work with the connectivity methods
   /// of the VTKMeshWriter class. If not zero, we will insert data in blocks of size blockStride.
   /// As first entry of each block we automatically insert (blockstride-1), which then describes
   /// the remaining number of entries in the block. For connectivity info these would be the node
   /// indices for the current element.
   template < typename entry_t, uint_t blockStride = 0u >
   class StreamAccessBuffer
   {
    public:
      StreamAccessBuffer( typename adios2::Variable< entry_t >::Span& buffer, adios2::Dims localDims )
      : buffer_( buffer )
      {
         if ( localDims.size() == 1 )
         {
            capacity_ = localDims[0];
         }
         else if ( localDims.size() == 2 )
         {
            capacity_ = localDims[0] * localDims[1];
         }
         else
         {
            WALBERLA_ABORT( " localDims.size() = " << localDims.size() );
         }

         // WALBERLA_LOG_INFO_ON_ROOT( "**** capacity of StreamAccessBuffer is " << capacity_ );
         // WALBERLA_LOG_INFO_ON_ROOT( "**** blockStride is " << blockStride );
      };

      StreamAccessBuffer& operator<<( const entry_t& value )
      {
         if constexpr ( blockStride > 0 )
         {
            if ( position_ % blockStride == 0 )
            {
               WALBERLA_ASSERT( position_ < capacity_ );
               buffer_[position_++] = blockStride - 1u;
            }
         }
         WALBERLA_ASSERT( position_ < capacity_ );
         buffer_[position_++] = value;
         return *this;
      }

    private:
      typename adios2::Variable< entry_t >::Span& buffer_;
      size_t                                      position_{ 0 };
      size_t                                      capacity_{ 0 };
   };

 private:
   /// object that remembers the functions we should export
   FEFunctionRegistry feFunctionRegistry_;

   /// exporting data will only work for functions associated with the mesh
   /// described be this PrimitiveStorage object
   std::shared_ptr< PrimitiveStorage > storage_{ nullptr };

   /// main ADIOS2 interface object
   //
   // ATTENTION: Order is important here! The adios_ member must come
   //            before all maps, so that it is destroyed after the
   //            writers in the maps!
   adios2::ADIOS adios_;

   /// we create one writer per level for P1 type functions
   std::map< uint_t, std::unique_ptr< AdiosWriterForP1 > > p1Writers_;

   /// we create one writer per level for P2 type functions
   std::map< uint_t, std::unique_ptr< AdiosWriterForP2 > > p2Writers_;

   /// path to output files
   std::string filePath_;

   /// basename of output files
   std::string fileBaseName_;

   /// type of engine to be used for export
   ///
   /// We will use the BP format, but the most recent "BP5" format is unsupported by ParaView 5.11.1.
   /// Thus, we still use BP4
   const std::string engineType_{ "BP4" };
};

} // namespace hyteg
