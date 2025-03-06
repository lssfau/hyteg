/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
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

#include "core/math/KahanSummation.h"

namespace hyteg {
namespace volumedofspace {

template < typename ValueType >
void VolumeDoFFunction< ValueType >::allocateMemory()
{
   if ( this->storage_->hasGlobalCells() )
   {
      // 3D

      // Create a data handling instance that handles the initialization, serialization, and deserialization of data.
      const auto dofDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >();

      // Create a data ID for all cells.
      this->storage_->addCellData( cellInnerDataID_, dofDataHandling, "VolumeDoFMacroCellData" );

      // Create a data handling instance that handles the initialization, serialization, and deserialization of data.
      const auto dofDataHandlingGL = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >();

      // Create three data IDs for all faces.
      this->storage_->addCellData( cellGhostLayerDataIDs_[0], dofDataHandling, "VolumeDoFMacroCellGL0Data" );
      this->storage_->addCellData( cellGhostLayerDataIDs_[1], dofDataHandling, "VolumeDoFMacroCellGL1Data" );
      this->storage_->addCellData( cellGhostLayerDataIDs_[2], dofDataHandling, "VolumeDoFMacroCellGL2Data" );
      this->storage_->addCellData( cellGhostLayerDataIDs_[3], dofDataHandling, "VolumeDoFMacroCellGL3Data" );

      // Allocate the DoFs.
      for ( auto& it : this->storage_->getCells() )
      {
         const auto pid  = it.first;
         const auto cell = it.second;

         // Fetching the FunctionMemory instance from each macro-face.
         FunctionMemory< ValueType >* functionMemory = cell->getData( cellInnerDataID_ );

         for ( uint_t level = this->minLevel_; level <= this->maxLevel_; level++ )
         {
            // Allocating the specified number of scalars on each micro-element for the entire macro-primitive.
            const auto numMacroLocalScalars = numScalarsPerPrimitive_.at( pid ) * levelinfo::num_microcells_per_cell( level );
            functionMemory->addData( level, numMacroLocalScalars, ValueType( 0 ) );
         }

         // Allocating ghost-layer memory only where necessary.
         for ( const auto& [localFaceID, npid] : cell->getIndirectNeighborCellIDsOverFaces() )
         {
            WALBERLA_UNUSED( npid );
            FunctionMemory< ValueType >* functionGLMemory = cell->getData( cellGhostLayerDataIDs_[localFaceID] );

            for ( uint_t level = this->minLevel_; level <= this->maxLevel_; level++ )
            {
               // TODO: adapt size here for mesh refinement
               const auto numGLScalars = numScalarsPerPrimitive_.at( pid ) * levelinfo::num_microfaces_per_face( level );
               functionGLMemory->addData( level, numGLScalars, ValueType( 0 ) );
            }
         }
      }
   }
   else
   {
      // 2D

      // Create a data handling instance that handles the initialization, serialization, and deserialization of data.
      const auto dofDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >();

      // Create a data ID for all faces.
      this->storage_->addFaceData( faceInnerDataID_, dofDataHandling, "VolumeDoFMacroFaceData" );

      // Create a data handling instance that handles the initialization, serialization, and deserialization of data.
      const auto dofDataHandlingGL = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >();

      // Create three data IDs for all faces.
      PrimitiveDataID< FunctionMemory< ValueType >, Face > fgldid0;
      PrimitiveDataID< FunctionMemory< ValueType >, Face > fgldid1;
      PrimitiveDataID< FunctionMemory< ValueType >, Face > fgldid2;

      this->storage_->addFaceData( fgldid0, dofDataHandlingGL, "VolumeDoFMacroFaceGL0Data" );
      this->storage_->addFaceData( fgldid1, dofDataHandlingGL, "VolumeDoFMacroFaceGL1Data" );
      this->storage_->addFaceData( fgldid2, dofDataHandlingGL, "VolumeDoFMacroFaceGL2Data" );

      faceGhostLayerDataIDs_[0] = fgldid0;
      faceGhostLayerDataIDs_[1] = fgldid1;
      faceGhostLayerDataIDs_[2] = fgldid2;

      // Allocate the DoFs.
      for ( auto& it : this->storage_->getFaces() )
      {
         const auto pid  = it.first;
         const auto face = it.second;

         // Fetching the FunctionMemory instance from each macro-face.
         FunctionMemory< ValueType >* functionMemory = face->getData( faceInnerDataID_ );

         for ( uint_t level = this->minLevel_; level <= this->maxLevel_; level++ )
         {
            // Allocating the specified number of scalars on each micro-element for the entire macro-primitive.
            const auto numMacroLocalScalars = numScalarsPerPrimitive_.at( pid ) * levelinfo::num_microfaces_per_face( level );
            functionMemory->addData( level, numMacroLocalScalars, ValueType( 0 ) );
         }

         // Allocating ghost-layer memory only where necessary.
         // The GL-level corresponds to the refinement level of the GL elements.
         // It has nothing to do with the macro refinement. This means depending on the neighboring macro refinement in case
         // of AMR applications, GL access must be performed using the respective level.
         for ( const auto& [localEdgeID, npids] : face->getIndirectTopLevelNeighborFaceIDsOverEdges() )
         {
            // The neighboring macros are either
            // a) on the same level (then there is only 1 macro)
            // b) on a coarser level (then there is also only 1 macro)
            // c) on a finer level (then there are multiple macros)
            // Since we only allocate a single ghost layer for each side regardless of the number of neighboring macros,
            // we can simply use any PID of the list and query the refinement level.
            auto npid = npids[0];

            FunctionMemory< ValueType >* functionGLMemory = face->getData( faceGhostLayerDataIDs_[localEdgeID] );

            // We potentially need to allocate more ghost-layers than minLevel and maxLevel indicate due to AMR refinement.
            // If the neighboring primitive is refined we must allocate from minLevel + 1 to maxLevel + 1.
            // If the neighboring primitive is coaser we must allocate from minLevel - 1 to maxLevel - 1.
            const auto primitiveLevel   = this->storage_->getRefinementLevel( pid );
            const auto glPrimitiveLevel = this->storage_->getRefinementLevel( npid );

            if ( glPrimitiveLevel == primitiveLevel - 1 || glPrimitiveLevel == primitiveLevel + 1 )
            {
               WALBERLA_CHECK_GREATER_EQUAL(
                   this->minLevel_,
                   1,
                   "AMR requires a minimum (micro-)refinement level of 1. Please change the minLevel parameter during function allocation." );
               WALBERLA_CHECK_GREATER_EQUAL(
                   this->maxLevel_,
                   1,
                   "AMR requires a minimum (micro-)refinement level of 1. Please change the maxLevel parameter during function allocation." );
            }
            else
            {
               WALBERLA_CHECK_EQUAL( glPrimitiveLevel,
                                     primitiveLevel,
                                     "2:1 balance does not seem to hold. Cannot work with that in VolumeDoFFunction." );
            }

            const uint_t allocationLevelOffset = glPrimitiveLevel - primitiveLevel;

            for ( uint_t level = this->minLevel_ + allocationLevelOffset; level <= this->maxLevel_ + allocationLevelOffset;
                  level++ )
            {
               // Handling ghost-layer size for AMR

               const auto numGLScalars = numScalarsPerPrimitive_.at( pid ) * levelinfo::num_microedges_per_edge( level );

               functionGLMemory->addData( level, numGLScalars, ValueType( 0 ) );
            }
         }
      }
   }
}

template < typename ValueType >
void VolumeDoFFunction< ValueType >::communicateNumScalarsPerPrimitive()
{
   walberla::mpi::BufferSystem bs( walberla::mpi::MPIManager::instance()->comm() );

   const auto neighboringRanks = this->storage_->getNeighboringVolumeRanksOfAllVolumes();

   bs.setReceiverInfo( neighboringRanks, true );

   // Complexity is not optimal - sending more data than necessary, but usually we do not need to send that often.
   for ( auto r : neighboringRanks )
   {
      bs.sendBuffer( r ) << numScalarsPerPrimitive_;
   }

   bs.sendAll();

   for ( auto i = bs.begin(); i != bs.end(); ++i )
   {
      std::map< PrimitiveID, uint_t > nScalarsNeighbor;
      i.buffer() >> nScalarsNeighbor;
      numScalarsPerPrimitive_.insert( nScalarsNeighbor.begin(), nScalarsNeighbor.end() );
   }
}

template < typename ValueType >
void VolumeDoFFunction< ValueType >::addPackInfos()
{
   for ( uint_t level = this->minLevel_; level <= this->maxLevel_; level++ )
   {
      communicators_[level]->addPackInfo( std::make_shared< VolumeDoFPackInfo< ValueType > >( this->getStorage(),
                                                                                              level,
                                                                                              numScalarsPerPrimitive_,
                                                                                              faceInnerDataID_,
                                                                                              cellInnerDataID_,
                                                                                              faceGhostLayerDataIDs_,
                                                                                              cellGhostLayerDataIDs_,
                                                                                              memoryLayout_ ) );
   }
}

template < typename ValueType >
VolumeDoFFunction< ValueType >::VolumeDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   const std::map< PrimitiveID, uint_t >&     numScalarsPerPrimitive,
                                                   indexing::VolumeDoFMemoryLayout            memoryLayout )

: Function< VolumeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, numScalarsPerPrimitive_( numScalarsPerPrimitive )
, memoryLayout_( memoryLayout )
, referenceCounter_( new internal::ReferenceCounter() )
{
   if ( this->storage_->getAdditionalHaloDepth() > 0 )
   {
      communicateNumScalarsPerPrimitive();
   }
   allocateMemory();
   addPackInfos();

   referenceCounter_->increaseRefs();
}

template < typename ValueType >
VolumeDoFFunction< ValueType >::VolumeDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   uint_t                                     numScalars,
                                                   indexing::VolumeDoFMemoryLayout            memoryLayout )

: Function< VolumeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, memoryLayout_( memoryLayout )
, referenceCounter_( new internal::ReferenceCounter() )
{
   numScalarsPerPrimitive_.clear();
   for ( auto pid : storage->getPrimitiveIDs() )
   {
      numScalarsPerPrimitive_[pid] = numScalars;
   }
   if ( this->storage_->getAdditionalHaloDepth() > 0 )
   {
      communicateNumScalarsPerPrimitive();
   }
   allocateMemory();
   addPackInfos();

   referenceCounter_->increaseRefs();
}

template < typename ValueType >
VolumeDoFFunction< ValueType >::~VolumeDoFFunction()
{
   referenceCounter_->decreaseRefs();
   if ( referenceCounter_->refs() <= 0 )
   {
      // There are no copies of this handle left. We can delete the allocated DoFs.
      deleteFunctionMemory();
   }
}

template < typename ValueType >
VolumeDoFFunction< ValueType >::VolumeDoFFunction( const VolumeDoFFunction< ValueType >& other )
: Function< VolumeDoFFunction< ValueType > >( other )
, numScalarsPerPrimitive_( other.numScalarsPerPrimitive_ )
, memoryLayout_( other.memoryLayout_ )
, faceInnerDataID_( other.faceInnerDataID_ )
, cellInnerDataID_( other.cellInnerDataID_ )
, faceGhostLayerDataIDs_( other.faceGhostLayerDataIDs_ )
, cellGhostLayerDataIDs_( other.cellGhostLayerDataIDs_ )
, referenceCounter_( other.referenceCounter_ )
{
   referenceCounter_->increaseRefs();
}

template < typename ValueType >
VolumeDoFFunction< ValueType >& VolumeDoFFunction< ValueType >::operator=( const VolumeDoFFunction< ValueType >& other )
{
   if ( this == &other )
   {
      return *this;
   }
   else if ( other.referenceCounter_ == referenceCounter_ )
   {
      if ( this->storage_->hasGlobalCells() )
      {
         WALBERLA_CHECK_EQUAL( cellInnerDataID_, other.cellInnerDataID_ )
      }
      else
      {
         WALBERLA_CHECK_EQUAL( faceInnerDataID_, other.faceInnerDataID_ )
      }
   }
   else
   {
      if ( this->storage_->hasGlobalCells() )
      {
         WALBERLA_CHECK_UNEQUAL( cellInnerDataID_, other.cellInnerDataID_ )
      }
      else
      {
         WALBERLA_CHECK_UNEQUAL( faceInnerDataID_, other.faceInnerDataID_ )
      }

      referenceCounter_->decreaseRefs();

      if ( referenceCounter_->refs() == 0 )
      {
         // There are no copies of this handle left. We can delete the allocated DoFs.
         deleteFunctionMemory();
      }

      numScalarsPerPrimitive_ = other.numScalarsPerPrimitive_;
      memoryLayout_           = other.memoryLayout_;
      faceInnerDataID_        = other.faceInnerDataID_;
      faceGhostLayerDataIDs_  = other.faceGhostLayerDataIDs_;
      cellInnerDataID_        = other.cellInnerDataID_;
      cellGhostLayerDataIDs_  = other.cellGhostLayerDataIDs_;
      referenceCounter_       = other.referenceCounter_;
      referenceCounter_->increaseRefs();
   }
   return *this;
}

template < typename ValueType >
void VolumeDoFFunction< ValueType >::communicate( uint_t level )
{
   if ( !this->storage_->hasGlobalCells() )
   {
      this->communicators_[level]->template startCommunication< Face, Face >();
      this->communicators_[level]->template endCommunication< Face, Face >();
   }
   else
   {
      this->communicators_[level]->template startCommunication< Cell, Cell >();
      this->communicators_[level]->template endCommunication< Cell, Cell >();
   }
}

/// \brief Assigns a linear combination of multiple VolumeDoFFunctions to this.
template < typename ValueType >
void VolumeDoFFunction< ValueType >::assign(
    const std::vector< ValueType >&                                                      scalars,
    const std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > >& functions,
    uint_t                                                                               level ) const
{
   WALBERLA_CHECK_EQUAL( scalars.size(),
                         functions.size(),
                         "VolumeDoFFunction< ValueType >::assign(): must pass same number of scalars and functions." )

   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& cellIt : this->getStorage()->getCells() )
      {
         const auto cellId = cellIt.first;
         const auto cell   = *cellIt.second;

         std::vector< ValueType* >                      srcPtrs( functions.size() );
         std::vector< indexing::VolumeDoFMemoryLayout > srcLayouts( functions.size() );
         for ( uint_t i = 0; i < functions.size(); i++ )
         {
            const auto f  = functions.at( i );
            srcPtrs[i]    = f.get().dofMemory( cellId, level );
            srcLayouts[i] = f.get().memoryLayout();
         }

         auto dstMem    = dofMemory( cellId, level );
         auto dstLayout = memoryLayout_;
         auto numDofs   = this->numScalarsPerPrimitive_.at( cellId );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( auto elementIdx : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  ValueType sum = 0;
                  for ( uint_t i = 0; i < functions.size(); i++ )
                  {
                     const auto s = scalars.at( i );

                     sum +=
                         s * srcPtrs[i][indexing::index(
                                 elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, srcLayouts[i] )];
                  }
                  dstMem[indexing::index(
                      elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, dstLayout )] = sum;
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         std::vector< ValueType* >                      srcPtrs( functions.size() );
         std::vector< indexing::VolumeDoFMemoryLayout > srcLayouts( functions.size() );
         for ( uint_t i = 0; i < functions.size(); i++ )
         {
            const auto f  = functions.at( i );
            srcPtrs[i]    = f.get().dofMemory( faceId, level );
            srcLayouts[i] = f.get().memoryLayout();
         }

         auto dstMem    = dofMemory( faceId, level );
         auto dstLayout = memoryLayout_;
         auto numDofs   = this->numScalarsPerPrimitive_.at( faceId );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  ValueType sum = 0;
                  for ( uint_t i = 0; i < functions.size(); i++ )
                  {
                     const auto s = scalars.at( i );

                     sum += s * srcPtrs[i][indexing::index(
                                    elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, srcLayouts[i] )];
                  }
                  dstMem[indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, dstLayout )] = sum;
               }
            }
         }
      }
   }
}

/// \brief Adds a scalar to this VolumeDoFFunction.
template < typename ValueType >
void VolumeDoFFunction< ValueType >::add( const ValueType scalar, uint_t level, DoFType flag ) const
{
   WALBERLA_UNUSED( flag );
   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& cellIt : this->getStorage()->getCells() )
      {
         const auto cellId = cellIt.first;
         const auto cell   = *cellIt.second;

         const auto mem     = dofMemory( cellId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( cellId );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( auto elementIdx : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx =
                      indexing::index( elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, layout );
                  mem[idx] += scalar;
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         auto       mem     = dofMemory( faceId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( faceId );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx = indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, layout );
                  mem[idx] += scalar;
               }
            }
         }
      }
   }
}

/// \brief Zero all function DoFs including those in ghost-layers
template < typename ValueType >
void VolumeDoFFunction< ValueType >::setToZero( uint_t level ) const
{
   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& cellIt : this->getStorage()->getCells() )
      {
         const auto cellID = cellIt.first;
         WALBERLA_CHECK( this->storage_->cellExistsLocally( cellID ),
                         "Cannot read/write DoF since macro-cell does not exists (locally)." );
         const Cell* cell = this->storage_->getCell( cellID );

         // zero cell DoFs
         FunctionMemory< ValueType >* fmem = cell->getData( cellInnerDataID_ );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         fmem->setToZero( level );

         // zero DoFs in ghost-layers
         for ( const auto& [localFaceID, npid] : cell->getIndirectNeighborCellIDsOverFaces() )
         {
            WALBERLA_UNUSED( npid );
            FunctionMemory< ValueType >* fmemGL = cell->getData( cellGhostLayerDataIDs_.at( localFaceID ) );
            fmemGL->setToZero( level );
         }
      }
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceID = faceIt.first;
         WALBERLA_CHECK( this->storage_->faceExistsLocally( faceID ),
                         "Cannot read/write DoF since macro-face does not exists (locally)." );
         const Face* face = this->storage_->getFace( faceID );

         // zero cell DoFs
         FunctionMemory< ValueType >* fmem = face->getData( faceInnerDataID_ );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         fmem->setToZero( level );

         // zero DoFs in ghost-layers
         for ( const auto& [localEdgeID, npid] : face->getIndirectTopLevelNeighborFaceIDsOverEdges() )
         {
            WALBERLA_UNUSED( npid );
            FunctionMemory< ValueType >* fmemGL = face->getData( faceGhostLayerDataIDs_.at( localEdgeID ) );
            fmemGL->setToZero( level );
         }
      }
   }
}

/// \brief Adds a linear combination of multiple VolumeDoFFunctions to this.
template < typename ValueType >
void VolumeDoFFunction< ValueType >::add(
    const std::vector< ValueType >&                                                      scalars,
    const std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > >& functions,
    uint_t                                                                               level ) const
{
   WALBERLA_CHECK_EQUAL( scalars.size(),
                         functions.size(),
                         "VolumeDoFFunction< ValueType >::add(): must pass same number of scalars and functions." )

   if ( this->storage_->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented" );
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         std::vector< ValueType* >                      srcPtrs( functions.size() );
         std::vector< indexing::VolumeDoFMemoryLayout > srcLayouts( functions.size() );
         for ( uint_t i = 0; i < functions.size(); i++ )
         {
            const auto f  = functions.at( i );
            srcPtrs[i]    = f.get().dofMemory( faceId, level );
            srcLayouts[i] = f.get().memoryLayout();
         }

         auto dstMem    = dofMemory( faceId, level );
         auto dstLayout = memoryLayout_;
         auto numDofs   = this->numScalarsPerPrimitive_.at( faceId );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  ValueType sum = 0;
                  for ( uint_t i = 0; i < functions.size(); i++ )
                  {
                     const auto s = scalars.at( i );

                     sum += s * srcPtrs[i][indexing::index(
                                    elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, srcLayouts[i] )];
                  }
                  dstMem[indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, dstLayout )] += sum;
               }
            }
         }
      }
   }
}

/// \brief Evaluates the dot product on all local DoFs. No communication is involved and the results may be different on each
/// process.
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::dotLocal( const VolumeDoFFunction< ValueType >& rhs, uint_t level ) const
{
   ValueType sum = 0;

   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& cellIt : this->getStorage()->getCells() )
      {
         const auto cellId = cellIt.first;
         const auto cell   = *cellIt.second;

         const auto mem     = dofMemory( cellId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( cellId );

         const auto otherMem    = rhs.dofMemory( cellId, level );
         const auto otherLayout = rhs.memoryLayout();

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( auto elementIdx : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx =
                      indexing::index( elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, layout );
                  const auto otherIdx = indexing::index(
                      elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, otherLayout );

                  sum += mem[idx] * otherMem[otherIdx];
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         const auto mem     = dofMemory( faceId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( faceId );

         const auto otherMem    = rhs.dofMemory( faceId, level );
         const auto otherLayout = rhs.memoryLayout();

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx = indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, layout );
                  const auto otherIdx =
                      indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, otherLayout );

                  sum += mem[idx] * otherMem[otherIdx];
               }
            }
         }
      }
   }

   return sum;
}

/// \brief swaps the content of one volumeDoFFunction with another.
template < typename ValueType >
void VolumeDoFFunction< ValueType >::swap( const VolumeDoFFunction< ValueType >& rhs, uint_t level ) const
{
   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& cellIt : this->getStorage()->getCells() )
      {
         const auto cellId = cellIt.first;
         const auto cell   = *cellIt.second;

         const auto mem     = dofMemory( cellId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( cellId );

         const auto otherMem    = rhs.dofMemory( cellId, level );
         const auto otherLayout = rhs.memoryLayout();

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( auto elementIdx : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx =
                      indexing::index( elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, layout );
                  const auto otherIdx = indexing::index(
                      elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, otherLayout );

                  std::swap( mem[idx], otherMem[otherIdx] );
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         const auto mem     = dofMemory( faceId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( faceId );

         const auto otherMem    = rhs.dofMemory( faceId, level );
         const auto otherLayout = rhs.memoryLayout();

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx = indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, layout );
                  const auto otherIdx =
                      indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, otherLayout );

                  std::swap( mem[idx], otherMem[otherIdx] );
               }
            }
         }
      }
   }
}
/// \brief Evaluates the (global) dot product. Involves communication and has to be called collectively.
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::dotGlobal( const VolumeDoFFunction< ValueType >& rhs, uint_t level ) const
{
   const auto dLocal = dotLocal( rhs, level );
   return walberla::mpi::allReduce( dLocal, walberla::mpi::SUM );
}

/// \brief Evaluates the sum of all local DoFs. No communication is involved and the results may be different on each
/// process.
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::sumLocal( uint_t level, bool absolute ) const
{
   walberla::math::KahanAccumulator< ValueType > sum;

   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& cellIt : this->getStorage()->getCells() )
      {
         const auto cellId = cellIt.first;
         const auto cell   = *cellIt.second;

         const auto mem     = dofMemory( cellId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( cellId );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( auto elementIdx : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx =
                      indexing::index( elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, layout );
                  sum += mem[idx];
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         const auto mem     = dofMemory( faceId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( faceId );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               if ( absolute )
               {
                  for ( uint_t dof = 0; dof < numDofs; dof++ )
                  {
                     const auto idx = indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, layout );
                     sum += std::abs( mem[idx] );
                  }
               }
               else
               {
                  for ( uint_t dof = 0; dof < numDofs; dof++ )
                  {
                     const auto idx = indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, layout );
                     sum += mem[idx];
                  }
               }
            }
         }
      }
   }

   return sum.get();
}

/// \brief Evaluates the (global) sum. Involves communication and has to be called collectively.
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::sumGlobal( uint_t level, bool absolute ) const
{
   const ValueType sLocal = sumLocal( level, absolute );
   return walberla::mpi::allReduce( sLocal, walberla::mpi::SUM );
}

// A generic local reduce operation to be used by getMaxDoFValue() and the like
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::reduceLocal( uint_t                                                    level,
                                                       const std::function< ValueType( ValueType, ValueType ) >& reduceOperation,
                                                       ValueType initialValue ) const
{
   ValueType result{ initialValue };

   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& cellIt : this->getStorage()->getCells() )
      {
         const auto cellId = cellIt.first;
         const auto cell   = *cellIt.second;

         const auto mem     = dofMemory( cellId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( cellId );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( auto elementIdx : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx =
                      indexing::index( elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, dof, numDofs, level, layout );
                  result = reduceOperation( result, mem[idx] );
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         const auto mem     = dofMemory( faceId, level );
         const auto layout  = memoryLayout_;
         const auto numDofs = this->numScalarsPerPrimitive_.at( faceId );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t dof = 0; dof < numDofs; dof++ )
               {
                  const auto idx = indexing::index( elementIdx.x(), elementIdx.y(), faceType, dof, numDofs, level, layout );
                  result         = reduceOperation( result, mem[idx] );
               }
            }
         }
      }
   }

   return result;
}

// Return the maximal value of the degrees of freedom of the function
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::getMaxDoFValue( uint_t level, bool mpiReduce ) const
{
   ValueType localMax = reduceLocal(
       level,
       []( ValueType oldMax, ValueType newCandidate ) { return std::max( oldMax, newCandidate ); },
       std::numeric_limits< ValueType >::lowest() );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

// Return the minimal value of the degrees of freedom of the function
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::getMinDoFValue( uint_t level, bool mpiReduce ) const
{
   ValueType localMin = reduceLocal(
       level,
       []( ValueType oldMin, ValueType newCandidate ) { return std::min( oldMin, newCandidate ); },
       std::numeric_limits< ValueType >::max() );

   ValueType globalMin = localMin;
   if ( mpiReduce )
   {
      globalMin = walberla::mpi::allReduce( localMin, walberla::mpi::MIN );
   }

   return globalMin;
}

// Return the maximal magnitude of the degrees of freedom of the function
template < typename ValueType >
ValueType VolumeDoFFunction< ValueType >::getMaxDoFMagnitude( uint_t level, bool mpiReduce ) const
{
   ValueType localMax = reduceLocal(
       level,
       []( ValueType oldMax, ValueType newCandidate ) { return std::max( oldMax, std::abs( newCandidate ) ); },
       ValueType( 0 ) );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
void VolumeDoFFunction< ValueType >::toVector( const VolumeDoFFunction< idx_t >&     numerator,
                                               const std::shared_ptr< VectorProxy >& vec,
                                               uint_t                                level,
                                               DoFType                               flag ) const
{
   if ( this->getStorage()->hasGlobalCells() )
   {
      // 3D
      for ( const auto& it : this->getStorage()->getCells() )
      {
         const auto cellID = it.first;
         const auto cell   = *it.second;

         const auto indices   = numerator.dofMemory( cellID, level );
         const auto dofs      = this->dofMemory( cellID, level );
         const auto memLayout = this->memoryLayout();
         const auto numDofs   = this->numScalarsPerPrimitive_.at( cellID );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const auto vectorIdx = indices[volumedofspace::indexing::index(
                      idxIt.x(), idxIt.y(), idxIt.z(), cellType, i, numDofs, level, memLayout )];
                  const auto value     = dofs[volumedofspace::indexing::index(
                      idxIt.x(), idxIt.y(), idxIt.z(), cellType, i, numDofs, level, memLayout )];
                  vec->setValue( uint_c( vectorIdx ), real_c( value ) );
               }
            }
         }
      }
   }
   else
   {
      // 2D
      for ( const auto& it : this->getStorage()->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         const auto indices   = numerator.dofMemory( faceID, level );
         const auto dofs      = this->dofMemory( faceID, level );
         const auto memLayout = this->memoryLayout();
         const auto numDofs   = this->numScalarsPerPrimitive_.at( faceID );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const auto vectorIdx =
                      indices[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )];
                  const auto value =
                      dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )];
                  vec->setValue( uint_c( vectorIdx ), real_c( value ) );
               }
            }
         }
      }
   }

   WALBERLA_UNUSED( flag );
}

template < typename ValueType >
void VolumeDoFFunction< ValueType >::fromVector( const VolumeDoFFunction< idx_t >&     numerator,
                                                 const std::shared_ptr< VectorProxy >& vec,
                                                 uint_t                                level,
                                                 DoFType                               flag ) const
{
   if ( this->getStorage()->hasGlobalCells() )
   {
      // 3D
      for ( const auto& it : this->getStorage()->getCells() )
      {
         const auto cellID = it.first;
         const auto cell   = *it.second;

         const auto indices   = numerator.dofMemory( cellID, level );
         auto       dofs      = this->dofMemory( cellID, level );
         const auto memLayout = this->memoryLayout();
         const auto numDofs   = this->numScalarsPerPrimitive_.at( cellID );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const auto vectorIdx = indices[volumedofspace::indexing::index(
                      idxIt.x(), idxIt.y(), idxIt.z(), cellType, i, numDofs, level, memLayout )];
                  const auto value     = vec->getValue( uint_c( vectorIdx ) );
                  dofs[volumedofspace::indexing::index(
                      idxIt.x(), idxIt.y(), idxIt.z(), cellType, i, numDofs, level, memLayout )] = ValueType( value );
               }
            }
         }
      }
   }
   else
   {
      // 2D
      for ( const auto& it : this->getStorage()->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         const auto indices   = numerator.dofMemory( faceID, level );
         auto       dofs      = this->dofMemory( faceID, level );
         const auto memLayout = this->memoryLayout();
         const auto numDofs   = this->numScalarsPerPrimitive_.at( faceID );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const auto vectorIdx =
                      indices[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )];
                  const auto value = vec->getValue( uint_c( vectorIdx ) );
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )] =
                      ValueType( value );
               }
            }
         }
      }
   }

   WALBERLA_UNUSED( flag );
}

// explicit instantiation
template class VolumeDoFFunction< double >;
template class VolumeDoFFunction< float >;
template class VolumeDoFFunction< int32_t >;
template class VolumeDoFFunction< int64_t >;

} // namespace volumedofspace

} // namespace hyteg
