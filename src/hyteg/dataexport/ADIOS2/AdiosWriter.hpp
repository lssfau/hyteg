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

class AdiosWriter : public FEFunctionWriter
{
 public:
   /// \param filePath          Path to directory where the files are stored
   /// \param fileBaseName      Basename of the vtk files
   /// \param storage           PrimitiveStorage associated with functions to export
   /// \param writeFrequency    Specifies the "frequency" of the exports see write()
   AdiosWriter( std::string                                filePath,
                std::string                                fileBaseName,
                const std::shared_ptr< PrimitiveStorage >& storage,
                const uint_t&                              writeFrequency = 1 )
   : AdiosWriter( filePath, fileBaseName, "", storage, writeFrequency )
   {}

   /// \param filePath          Path to directory where the files are stored
   /// \param fileBaseName      Basename of the vtk files
   /// \param configFile        Name of a file in XML or YAML format with runtime configuration parameters for ADIOS2
   /// \param storage           PrimitiveStorage associated with functions to export
   /// \param writeFrequency    Specifies the "frequency" of the exports see write()
   AdiosWriter( std::string                                filePath,
                std::string                                fileBaseName,
                std::string                                configFile,
                const std::shared_ptr< PrimitiveStorage >& storage,
                const uint_t&                              writeFrequency = 1 )
   : storage_( storage )
   , filePath_( filePath )
   , fileBaseName_( fileBaseName )
   , configFile_( configFile )
   {}

   /// The destructor takes care of ... nothing
   ~AdiosWriter() {
      WALBERLA_LOG_INFO_ON_ROOT( "D'tor of AdiosWriter called" );
   }

   /// Add an FE Function to became part of the next dataexport phase
   template < template < typename > class func_t, typename value_t >
   inline void add( const func_t< value_t >& function )
   {
      if constexpr ( !std::is_same_v< func_t< value_t >, P1Function< value_t > > &&
                     !std::is_same_v< func_t< value_t >, P2Function< value_t > > &&
                     !std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > &&
                     !std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
      {
         WALBERLA_ABORT( "AdiosWriter only supports P1 and P2 type functions!" );
      }

      feFunctionRegistry_.add( function );
   }

   /// Set parameter specified by string key to value specified by string value
   ///
   /// Currently not supported for AdiosWriter
   void setParameter( const std::string& key, const std::string& value ) override final
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter::setParameter() does not perform any action!" );
   }

   void write( const uint_t level, const uint_t timestep = 0 ) override final;

   /// Class that wraps an Adios span such that we can insert data with operator<<
   ///
   /// The blockStride template argument is designed to work with the connectivity methods
   /// of the VTKMeshWriter class. If not zero, we will insert data in blocks of size blockStride.
   /// As first entry of each block we automatically insert (blockstride-1), which then described
   /// the remaining number of entries in the block. For connectivity info these would be the node
   /// indices for the current element.
   template < typename entry_t, uint_t blockStride = 0 >
   class StreamAccessBuffer
   {
    public:
      StreamAccessBuffer( typename adios2::Variable< entry_t >::Span& buffer, adios2::Dims localDims )
      : buffer_( buffer )
      {
         WALBERLA_ASSERT( localDims.size() == 2 );
         capacity_ = localDims[0] * localDims[1];
         WALBERLA_LOG_INFO_ON_ROOT( "**** capacity of StreamAccessBuffer is " << capacity_ );
         WALBERLA_LOG_INFO_ON_ROOT( "**** blockStride is " << blockStride );
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

   mutable std::map< uint_t, std::unique_ptr< AdiosWriterForP1 > > p1Writers_;
   // std::map< uint_t, AdiosWriterForP2 > p2Writers_;

   std::string filePath_;
   std::string fileBaseName_;
   std::string configFile_;
};

} // namespace hyteg
