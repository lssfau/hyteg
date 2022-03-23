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

#pragma once

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"
#include "hyteg/volumedofspace/VolumeDoFPackInfo.hpp"

namespace hyteg {
namespace volumedofspace {

/// \brief Parallel data structure for degrees of freedom that live in volumes.
///
/// The VolumeDoFFunction has the following properties:
///  - a specified number of scalars is allocated per micro-volume
///  - in 2D DoFs are allocated on macro-faces only (micro-volumes are micro-faces)
///  - in 3D DoFs are allocated on macro-cells only (micro-volumes are micro-cells)
///  - ghost-layers are allocated so that each inner micro-volume of a macro-primitives has neighboring DoFs interfaced over
///    all facets (dim - 1, that is edges in 2D (3 facets) and faces in 3D (4 facets))
///  - the number of scalars allocated on each macro-volume can be specified individually to allow for the implementation of e.g.
///    adaptive p-refinement, the scalars on the ghost-layers are adapted accordingly
///  - there is no "position" assigned to the DoFs - meaning to the scalars has to be given by external data structures
template < typename ValueType >
class VolumeDoFFunction : public Function< VolumeDoFFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = VolumeDoFFunction< VType >;

   /// \brief Allocates a VolumeDoFFunction.
   ///
   /// \param name                   name of this function
   /// \param storage                PrimitiveStorage object
   /// \param minLevel               smallest h-refinement level to allocate
   /// \param maxLevel               largest h-refinement level to allocate
   /// \param numScalarsPerPrimitive map of (local) PrimitiveIDs to the number of DoFs per micro-element that shall be allocated
   /// \param memoryLayout           Structure-of-Arrays or Array-of-Structures
   VolumeDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      const std::map< PrimitiveID, uint_t >&     numScalarsPerPrimitive,
                      indexing::VolumeDoFMemoryLayout            memoryLayout );

   /// \brief Allocates a VolumeDoFFunction.
   ///
   /// \param name         name of this function
   /// \param storage      PrimitiveStorage object
   /// \param minLevel     smallest h-refinement level to allocate
   /// \param maxLevel     largest h-refinement level to allocate
   /// \param numScalars   number of scalars allocated per element (same for all macro-volumes)
   /// \param memoryLayout Structure-of-Arrays or Array-of-Structures
   VolumeDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      uint_t                                     numScalars,
                      indexing::VolumeDoFMemoryLayout            memoryLayout );

   /// \brief Updates ghost-layers.
   void communicate( uint_t level );

   /// \brief Assigns a linear combination of multiple VolumeDoFFunctions to this.
   void assign( const std::vector< ValueType >&                                                      scalars,
                const std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > >& functions,
                uint_t                                                                               level );

   /// \brief Evaluates the dot product on all local DoFs. No communication is involved and the results may be different on each
   /// process.
   ValueType dotLocal( const VolumeDoFFunction< ValueType >& rhs, uint_t level ) const;

   /// \brief Evaluates the (global) dot product. Involves communication and has to be called collectively.
   ValueType dotGlobal( const VolumeDoFFunction< ValueType >& rhs, uint_t level ) const;

   /// \brief Returns a pointer to the array that stores all degrees of freedom.
   ///
   /// You better know what you are doing when calling this ...
   ///
   /// \param primitiveID macro-primitive (volume primitive) to get the data of
   /// \return pointer to the array that stores the data
   ValueType* dofMemory( PrimitiveID primitiveID, uint_t level ) const
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_CHECK( storage_->cellExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-cell does not exists (locally)." );
         auto fmem = storage_->getCell( primitiveID )->template getData( cellInnerDataID_ );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         auto data = fmem->getPointer( level );
         return data;
      }
      else
      {
         WALBERLA_CHECK( storage_->faceExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-face does not exists (locally)." );
         auto fmem = storage_->getFace( primitiveID )->template getData( faceInnerDataID_ );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         auto data = fmem->getPointer( level );
         return data;
      }
   }

   /// \brief Returns a pointer to the array that stores the ghost-layer data.
   ///
   /// You better know what you are doing when calling this ...
   ///
   /// \param primitiveID macro-primitive (volume primitive) to get the data of
   /// \return pointer to the array that stores the data
   ValueType* glMemory( PrimitiveID primitiveID, uint_t level, uint_t glId ) const
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_CHECK( storage_->cellExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-cell does not exists (locally)." );
         FunctionMemory< ValueType >* fmem =
             storage_->getCell( primitiveID )->template getData( cellGhostLayerDataIDs_.at( glId ) );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         ValueType* data = fmem->getPointer( level );
         return data;
      }
      else
      {
         WALBERLA_CHECK( storage_->faceExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-face does not exists (locally)." );
         FunctionMemory< ValueType >* fmem =
             storage_->getFace( primitiveID )->template getData( faceGhostLayerDataIDs_.at( glId ) );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         ValueType* data = fmem->getPointer( level );
         return data;
      }
   }

   /// \brief Read access to a degree of freedom on a macro-face (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param faceType    face type of the micro-volume
   /// \param level       refinement level
   /// \return value of the DoF
   ValueType
       dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, facedof::FaceType faceType, uint_t level ) const
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), faceType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   /// \brief Access to a degree of freedom on a macro-face (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param faceType    face type of the micro-volume
   /// \param level       refinement level
   /// \return reference to the DoF
   ValueType& dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, facedof::FaceType faceType, uint_t level )
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), faceType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   /// \brief Read access to a degree of freedom on a macro-cell (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param cellType    cell type of the micro-volume
   /// \param level       refinement level
   /// \return value of the DoF
   ValueType
       dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, celldof::CellType cellType, uint_t level ) const
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), idx.z(), cellType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   /// \brief Access to a degree of freedom on a macro-cell (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param cellType    cell type of the micro-volume
   /// \param level       refinement level
   /// \return reference to the DoF
   ValueType& dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, celldof::CellType cellType, uint_t level )
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), idx.z(), cellType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   indexing::VolumeDoFMemoryLayout memoryLayout() const { return memoryLayout_; }

 private:
   using Function< VolumeDoFFunction< ValueType > >::communicators_;

   void allocateMemory();

   /// \brief Updates the map holding the number of scalars per primitive with neighborhood information.
   void communicateNumScalarsPerPrimitive();

   void addPackInfos();

   std::string name_;

   std::shared_ptr< PrimitiveStorage > storage_;

   uint_t minLevel_;
   uint_t maxLevel_;

   std::map< PrimitiveID, uint_t > numScalarsPerPrimitive_;

   indexing::VolumeDoFMemoryLayout memoryLayout_;

   /// Each data ID stores the inner scalar fields for all levels.
   PrimitiveDataID< FunctionMemory< ValueType >, Face > faceInnerDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell > cellInnerDataID_;

   /// One data ID per ghost-layer. Should be up to 3 in 2D and 4 in 3D.
   /// Maps from the local macro-edge ID (for 2D) or local macro-face ID (for 3D) to the respective ghost-layer memory.
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Face > > faceGhostLayerDataIDs_;
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Cell > > cellGhostLayerDataIDs_;
};

/// \brief Given an affine coordinate (in computational space) this function locates the micro-element in the macro-face and
/// the coordinates (in reference space) local to the micro-element.
///
/// \param level            [in] refinement level
/// \param face             [in] macro-face object
/// \param coordinates      [in] the point that the local element shall include
/// \param elementIndex     [out] logical index of the micro-element
/// \param faceType         [out] face type (up- or down-face / blue or grey face)
/// \param localCoordinates [out] coordinates of the queried point translated to element local coordinates (reference space)
template < typename ValueType >
inline void getLocalElementFromCoordinates( uint_t                  level,
                                            const Face&             face,
                                            const Point2D           coordinates,
                                            hyteg::indexing::Index& elementIndex,
                                            facedof::FaceType&      faceType,
                                            Point2D&                localCoordinates )
{
   // TODO: precompute and store trafo matrix from affine to ref space in macro-face
   //       and/or write method in macro-face that performs the trafo altogether

   // Transform absolute coordinates to macro element relative coordinates
   Matrix2r A;
   A( 0, 0 )          = ( face.getCoordinates()[1] - face.getCoordinates()[0] )[0];
   A( 0, 1 )          = ( face.getCoordinates()[2] - face.getCoordinates()[0] )[0];
   A( 1, 0 )          = ( face.getCoordinates()[1] - face.getCoordinates()[0] )[1];
   A( 1, 1 )          = ( face.getCoordinates()[2] - face.getCoordinates()[0] )[1];
   Matrix2r transform = A.adj();
   transform *= 1.0 / A.det();

   Point2D x( { coordinates[0] - face.getCoordinates()[0][0], coordinates[1] - face.getCoordinates()[0][1] } );

   Point2D xRelMacro = transform.mul( x );

   // Determine lower-left corner index of the quad where the evaluation point lies in
   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );
   real_t hInv    = walberla::real_c( rowsize - 1 );
   real_t h       = walberla::real_c( 1.0 / hInv );
   int    binX    = static_cast< int >( std::floor( xRelMacro[0] * ( rowsize - 1 ) ) );
   int    binY    = static_cast< int >( std::floor( xRelMacro[1] * ( rowsize - 1 ) ) );

   if ( binX < 0 )
   {
      binX = 0;
   }

   if ( binY < 0 )
   {
      binY = 0;
   }

   if ( binX >= static_cast< int >( rowsize - 1 ) )
   {
      binX = static_cast< int >( rowsize - 2 );
   }

   if ( binY >= static_cast< int >( rowsize - 1 ) )
   {
      binY = static_cast< int >( rowsize - 2 );
   }

   if ( binX + binY >= static_cast< int >( rowsize - 1 ) )
   {
      int binXDec = ( binX + binY - static_cast< int >( rowsize - 2 ) ) / 2;
      binXDec += ( binX + binY - static_cast< int >( rowsize - 2 ) ) % 2;
      int binYDec = ( binX + binY - static_cast< int >( rowsize - 2 ) ) / 2;
      binX -= binXDec;
      binY -= binYDec;
   }

   hyteg::indexing::Index index;
   index.x() = uint_c( binX );
   index.y() = uint_c( binY );

   WALBERLA_ASSERT_LESS( index.x(), rowsize - 1 );
   WALBERLA_ASSERT_LESS( index.y(), rowsize - 1 );
   WALBERLA_ASSERT_LESS( index.x() + index.y(), rowsize - 1, "index.x(): " << index.x() << ", index.y()" << index.y() );

   Point2D localCoords;
   localCoords[0] = xRelMacro[0] - index.x() * h;
   localCoords[1] = xRelMacro[1] - index.y() * h;
   localCoords *= hInv;

   // decide if up or down triangle
   // clamp to macro-face if the corresponding down-triangle would be out of the macro-face
   // otherwise check floating point distance
   bool upTriangle = ( index.x() + index.y() == rowsize - 2 ) || ( localCoordinates[0] + localCoordinates[1] <= 1.0 );

   elementIndex     = index;
   faceType         = upTriangle ? facedof::FaceType::GRAY : facedof::FaceType::BLUE;
   localCoordinates = localCoords;
}

} // namespace volumedofspace
} // namespace hyteg