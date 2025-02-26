/*
 * Copyright (c) 2023-2024 Marcus Mohr, Roman Freissler.
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

#include "hyteg/HytegDefinitions.hpp"

#ifdef HYTEG_BUILD_WITH_ADIOS2

#include <adios2.h>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP1.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP2.hpp"
#include "hyteg/dataexport/FEFunctionWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKMeshWriter.hpp"
#include "hyteg/functions/FEFunctionRegistry.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

using adiosHelpers::adiostype_t;

class AdiosWriterForP1;
class AdiosWriterForP2;

/// Class to export FEFunction data using the ADIOS2 library
///
/// This class allows to write FEFunction data to the filesystem using the ADIOS2 library.
/// Output is using the BP format (see engineType_ member for version). The class currently
/// only supports the following types of functions:
/// - P1Function and P1VectorFunction
/// - P2Function and P2VectorFunction
/// - P2P1TaylorHoodFunction
///
/// which must be using real_t as their value type.
///
/// \note For configuring writer behaviour through an XML or YAML configuration file, please
///       note that we need to create one adios2::IO object for each output level and for
///       each function family. Their names are given by
///       - AdiosWriterP1-level<level number>
///       - AdiosWriterP2-level<level number>
class AdiosWriter : public FEFunctionWriter< AdiosWriter >
{
 public:
#ifdef HYTEG_BUILD_WITH_MPI
   /// \param filePath          Path to directory where the BP files are stored
   /// \param fileBaseName      Basename for output BP file (which is acutally a folder)
   /// \param storage           PrimitiveStorage associated with functions to export
   /// \param comm              MPI Communicator, defaults to the HyTeG standard communicator
   AdiosWriter( const std::string&                         filePath,
                const std::string&                         fileBaseName,
                const std::shared_ptr< PrimitiveStorage >& storage,
                MPI_Comm                                   comm = walberla::MPIManager::instance()->comm() )
   : AdiosWriter( filePath, fileBaseName, "", storage, comm )
   {}

   /// \param filePath          Path to directory where the BP files are stored
   /// \param fileBaseName      Basename for output BP file (which is acutally a folder)
   /// \param configFile        Name of a file in XML or YAML format with runtime configuration parameters for ADIOS2
   /// \param storage           PrimitiveStorage associated with functions to export
   /// \param comm              MPI Communicator, defaults to the HyTeG standard communicator
   AdiosWriter( const std::string&                         filePath,
                const std::string&                         fileBaseName,
                const std::string&                         configFile,
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
   AdiosWriter( const std::string& filePath, const std::string& fileBaseName, const std::shared_ptr< PrimitiveStorage >& storage )
   : AdiosWriter( filePath, fileBaseName, "", storage )
   {}

   /// \param filePath          Path to directory where the BP files are stored
   /// \param fileBaseName      Basename for output BP file (which is acutally a folder)
   /// \param configFile        Name of a file in XML or YAML format with runtime configuration parameters for ADIOS2
   /// \param storage           PrimitiveStorage associated with functions to export
   AdiosWriter( const std::string&                         filePath,
                const std::string&                         fileBaseName,
                const std::string&                         configFile,
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

   /// Add an FE Function to become part of the next dataexport phase
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

      if ( firstWriteDidHappen_ )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter class does not support adding functions after the first write.\n"
                                       << "--> Ignoring function '" << function.getFunctionName() << "'!" );
      }
      else
      {
         feFunctionRegistry_.add( function );
      }
   }

   /// Add a single-value attribute to ADIOS2 output
   ///
   /// Attributes are added to a map and will be written at the first
   /// invocation of write() on the AdiosWriter object
   /// \param key               string for map
   /// \param attrValue         scalar variable with ADIOS2-compatible data type or bool
   template < typename attribute_t >
   inline void addAttribute( const std::string& key, attribute_t attrValue )
   {
      // Check if attribute has been collected already
      if ( userDefinedAttributes_.count( key ) > 0 )
      {
         WALBERLA_ABORT( "The attribute " << key << " has been defined multiple times" );
      }

      // Ensure that function was instantiated for a supported attribute type
      if constexpr ( !(adiosHelpers::isAdiosDataType< attribute_t, adiostype_t >::value ||
                       std::is_same_v< attribute_t, const char* >) )
      {
         WALBERLA_ABORT( "The provided attribute " << key << " does not hold a valid datatype for ADIOS2." );
      }

      if ( firstWriteDidHappen_ )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter class does not support adding attributes after the first write.\n"
                                       << "--> Ignoring (key,value) = ('" << key << "','" << attrValue << "')" );
      }
      else
      {
         if constexpr ( std::is_same_v< attribute_t, const char* > )
            userDefinedAttributes_[key] = std::string{ attrValue };
         else
            userDefinedAttributes_[key] = attrValue;
      }
   }

   /// Set all attributes to be added to ADIOS2 output
   ///
   /// Attributes will be written at the first
   /// invocation of write() on the AdiosWriter object
   /// \param key               map of attributes
   inline void setAllAttributes( std::map< std::string, adiostype_t > userDefinedAttributes )
   {
      userDefinedAttributes_ = userDefinedAttributes;
   }

   /// Set parameter specified by string key to value specified by string value
   ///
   /// \note For consistency reasons the method will refuse to change parameter values after the first
   ///       invocation of write() on the AdiosWriter object. In this case it will print a warning and
   ///       ignore the request.
   void setParameter( const std::string& key, const std::string& value )
   {
      if ( firstWriteDidHappen_ )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriter::setParameter() only works before the first write()!\n"
                                       << "--> Ignoring (key,value) = ('" << key << "','" << value << "')" );
      }
      else
      {
         userProvidedParameters_[key] = value;
      }
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
      };

      template < typename input_t >
      StreamAccessBuffer& operator<<( const input_t& value )
      {
         if constexpr ( blockStride > 0 )
         {
            if ( position_ % blockStride == 0 )
            {
               WALBERLA_ASSERT( position_ < capacity_ );
               buffer_[position_++] = static_cast< entry_t >( blockStride - 1u );
            }
         }
         WALBERLA_ASSERT( position_ < capacity_ );
         buffer_[position_++] = static_cast< entry_t >( value );
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
   /// described by this PrimitiveStorage object
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

   /// (key,value) pairs for IO parameters provided by user via setParameter() method
   std::map< std::string, std::string > userProvidedParameters_;

   /// (key,value) pairs for adding attributes to ADIOS2 output
   std::map< std::string, adiostype_t > userDefinedAttributes_;

   /// remember if we already had a write() episode
   bool firstWriteDidHappen_ = false;

   /// type of engine to be used for export
   ///
   /// Currently the most recent format is "BP5" which is supported by ParaView
   /// from at least 5.12.1
   inline static const std::string engineType_{ "BP5" };
};

} // namespace hyteg

#endif
