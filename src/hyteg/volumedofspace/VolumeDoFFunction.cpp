/*
* Copyright (c) 2017-2022 Nils Kohl.
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


#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {
namespace volumedofspace {

template < typename ValueType >
void VolumeDoFFunction< ValueType >::allocateMemory()
{
   if ( storage_->hasGlobalCells() )
   {
      // 3D
      WALBERLA_ABORT( "DG not implemented in 3D" );
   }
   else
   {
      // 2D

      // Create a data handling instance that handles the initialization, serialization, and deserialization of data.
      const auto dofDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >();

      // Create a data ID for all faces.
      storage_->template addFaceData( faceInnerDataID_, dofDataHandling, "VolumeDoFMacroFaceData" );

      // Allocate the DoFs.
      for ( auto& it : storage_->getFaces() )
      {
         const auto pid  = it.first;
         const auto face = it.second;

         // Fetching the FunctionMemory instance from each macro-face.
         FunctionMemory< ValueType >* functionMemory = face->template getData( faceInnerDataID_ );

         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            // Allocating the specified number of scalars on each micro-element for the entire macro-primitive.
            const auto numMacroLocalScalars = numScalarsPerPrimitive_.at( pid ) * levelinfo::num_microfaces_per_face( level );
            functionMemory->addData( level, numMacroLocalScalars, ValueType( 0 ) );
         }
      }
   }
}

template < typename ValueType >
VolumeDoFFunction< ValueType >::VolumeDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   const std::map< PrimitiveID, uint_t >&     numScalarsPerPrimitive,
                                                   VolumeDoFMemoryLayout                      memoryLayout )

: Function< VolumeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, name_( name )
, storage_( storage )
, minLevel_( minLevel )
, maxLevel_( maxLevel )
, numScalarsPerPrimitive_( numScalarsPerPrimitive )
, memoryLayout_( memoryLayout )
{
   allocateMemory();
}

template < typename ValueType >
VolumeDoFFunction< ValueType >::VolumeDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   uint_t                                     numScalars,
                                                   VolumeDoFMemoryLayout                      memoryLayout )

: Function< VolumeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, name_( name )
, storage_( storage )
, minLevel_( minLevel )
, maxLevel_( maxLevel )
, memoryLayout_( memoryLayout )
{
   numScalarsPerPrimitive_.clear();
   for ( auto pid : storage->getPrimitiveIDs() )
   {
      numScalarsPerPrimitive_[pid] = numScalars;
   }

   allocateMemory();
}

/// explicit instantiation
template class VolumeDoFFunction< double >;
template class VolumeDoFFunction< float >;
template class VolumeDoFFunction< int32_t >;
template class VolumeDoFFunction< int64_t >;

} // namespace volumedofspace

} // namespace hyteg
