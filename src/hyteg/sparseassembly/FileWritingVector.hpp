/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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
#include <fstream>
#include <iostream>
#include <map>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/sparseassembly/MapVector.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

template < class FunctionType >
class FileWritingVector : public VectorProxy
{
 public:
   using indexType = typename FunctionType::template FunctionType< idx_t >;

   FileWritingVector( const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     level,
                      const indexType&                           enumerator,
                      std::shared_ptr< FunctionType >            function   = nullptr,
                      uint_t                                     bufferSize = 1024 )
   : bufferSize_( bufferSize )
   {
      mapVector_ = std::make_shared< MapVector >();

      if ( function == nullptr )
      {
         function = std::make_shared< FunctionType >( "FileWritingVector fct", storage, level, level );
      }

      function->toVector( enumerator, mapVector_, level, All );

      localMapSize_  = mapVector_->getSize();
      globalMapSize_ = walberla::mpi::allReduce( localMapSize_, walberla::mpi::SUM );
      floatSize_     = sizeof( real_t );
   }
   virtual ~FileWritingVector() {}

   void setValue( uint_t idx, real_t value ) override { mapVector_->setValue( idx, value ); }

   real_t getValue( uint_t idx ) const override { return mapVector_->getValue( idx ); }

   void readFromFile( std::string filename )
   {
      WALBERLA_MPI_WORLD_BARRIER();

      uint_t newLocalMapSize  = mapVector_->getSize();
      uint_t newGlobalMapSize = walberla::mpi::allReduce( newLocalMapSize, walberla::mpi::SUM );

      if ( newGlobalMapSize != globalMapSize_ )
      {
         WALBERLA_ABORT( "Number of DoFs has changed. Did you add additional values after the constructor call?" );
      }

      WALBERLA_ROOT_SECTION()
      {
         std::ifstream file( filename, std::ios::binary | std::ios::in );

         if ( file.is_open() )
         {
            file.seekg( 0, std::ios::end );

            std::streampos size = file.tellg();

            if ( size < globalMapSize_ * floatSize_ )
            {
               WALBERLA_ABORT( "Input file too small!" );
            }
         }
         else
         {
            WALBERLA_ABORT( "Error opening file " << filename << "!" );
         }

         file.close();
      }

      WALBERLA_MPI_WORLD_BARRIER();

      // read from the file
      {
         std::ifstream file( filename, std::ios::binary | std::ios::in );

         if ( file.is_open() )
         {
            std::vector< real_t > readBuffer( bufferSize_, real_c( 0 ) );
            uint_t                currentIndex = 0;
            uint_t                initIndex    = 0;
            bool                  first        = true;
            auto                  valueMap     = mapVector_->getValueMap();
            for ( auto& [key, value] : ( *valueMap ) )
            {
               if ( first )
               {
                  initIndex    = key;
                  currentIndex = key;
                  first        = false;
               }
               else
               {
                  uint_t newDataSize = ( currentIndex - initIndex + 1 );
                  if ( key > currentIndex + 1 || newDataSize >= bufferSize_ )
                  {
                     file.seekg( static_cast< std::streamsize >( initIndex * floatSize_ ), std::ios::beg );
                     file.read( reinterpret_cast< char* >( readBuffer.data() ),
                                static_cast< std::streamsize >( floatSize_ * newDataSize ) );
                     for ( uint_t i = 0; i < newDataSize; i++ )
                     {
                        ( *valueMap )[initIndex + i] = readBuffer[i];
                     }

                     initIndex    = key;
                     currentIndex = key;
                  }
                  else
                  {
                     currentIndex = key;
                  }
               }
            }

            if ( !first )
            {
               uint_t finalDataSize = ( currentIndex - initIndex + 1 );
               file.seekg( static_cast< std::streamsize >( initIndex * floatSize_ ), std::ios::beg );
               file.read( reinterpret_cast< char* >( readBuffer.data() ),
                          static_cast< std::streamsize >( floatSize_ * finalDataSize ) );
               for ( uint_t i = 0; i < finalDataSize; i++ )
               {
                  ( *valueMap )[initIndex + i] = readBuffer[i];
               }
            }
         }
         else
         {
            WALBERLA_ABORT( "Error opening file " << filename << "!" );
         }

         file.close();
      }

      WALBERLA_MPI_WORLD_BARRIER();
   }

   void writeToFile( std::string filename )
   {
      WALBERLA_MPI_WORLD_BARRIER();

      uint_t newLocalMapSize  = mapVector_->getSize();
      uint_t newGlobalMapSize = walberla::mpi::allReduce( newLocalMapSize, walberla::mpi::SUM );

      if ( newGlobalMapSize != globalMapSize_ )
      {
         WALBERLA_ABORT( "Number of DoFs has changed. Did you add additional values after the constructor call?" );
      }

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file( filename, std::ios::binary | std::ios::out );

         if ( file.is_open() )
         {
            // Preallocate space
            std::vector< real_t > writeBuffer( bufferSize_, real_c( 0 ) );
            uint_t                numberBufferWrites = globalMapSize_ / bufferSize_;
            uint_t                remainingWriteSize = globalMapSize_ - numberBufferWrites * bufferSize_;

            file.seekp( 0, std::ios::beg );

            for ( uint_t k = 0; k < numberBufferWrites; k++ )
            {
               file.write( reinterpret_cast< char* >( writeBuffer.data() ),
                           static_cast< std::streamsize >( floatSize_ * writeBuffer.size() ) );
            }
            file.write( reinterpret_cast< char* >( writeBuffer.data() ),
                        static_cast< std::streamsize >( floatSize_ * remainingWriteSize ) );
         }
         else
         {
            WALBERLA_ABORT( "Error creating file " << filename << "!" );
         }

         file.close();
      }

      WALBERLA_MPI_WORLD_BARRIER();

      // write to the file
      {
         std::ofstream file( filename, std::ios::binary | std::ios::out | std::ios::in );

         if ( file.is_open() )
         {
            std::vector< real_t > writeBuffer;
            uint_t                currentIndex = 0;
            uint_t                initIndex    = 0;
            auto                  valueMap     = mapVector_->getValueMap();
            for ( auto& [key, value] : ( *valueMap ) )
            {
               if ( writeBuffer.size() == 0 )
               {
                  currentIndex = key;
                  initIndex    = key;
                  writeBuffer.push_back( value );
               }
               else
               {
                  if ( key > currentIndex + 1 || writeBuffer.size() > bufferSize_ )
                  {
                     file.seekp( static_cast< std::streamsize >( initIndex * floatSize_ ), std::ios::beg );
                     file.write( reinterpret_cast< char* >( writeBuffer.data() ),
                                 static_cast< std::streamsize >( floatSize_ * writeBuffer.size() ) );

                     writeBuffer.clear();
                     currentIndex = key;
                     initIndex    = key;
                     writeBuffer.push_back( value );
                  }
                  else
                  {
                     currentIndex = key;
                     writeBuffer.push_back( value );
                  }
               }
            }

            file.seekp( static_cast< std::streamsize >( initIndex * floatSize_ ), std::ios::beg );
            file.write( reinterpret_cast< char* >( writeBuffer.data() ),
                        static_cast< std::streamsize >( floatSize_ * writeBuffer.size() ) );
         }
         else
         {
            WALBERLA_ABORT( "Error opening file " << filename << "!" );
         }

         file.close();
      }

      WALBERLA_MPI_WORLD_BARRIER();
   }

 private:
   std::shared_ptr< MapVector > mapVector_;

   uint_t localMapSize_;
   uint_t globalMapSize_;
   uint_t floatSize_;
   uint_t bufferSize_;
};

} // namespace hyteg
