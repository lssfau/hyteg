/*
* Copyright (c) 2017-2021 Nils Kohl.
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

namespace hyteg {
namespace volumedofspace {

/// \brief Defines the memory layout for a VolumeDoFFunction.
///
/// structure-of-arrays (SoA): the innermost variable determines the index of the current micro-volume
/// array-of-structures (AoS): the innermost variable determines the DoF of the current micro-volume
enum class VolumeDoFMemoryLayout
{
   SoA,
   AoS
};

namespace indexing {

/// \brief Converts the refinement level to the 'width' (number of micro-edges on the macro-edge) of the volume.
inline constexpr uint_t levelToWidth( uint_t level )
{
   return levelinfo::num_microedges_per_edge( level );
}

/// \brief Returns the array-index of the specified 2D volume DoF.
///
/// \param x         logical x-coordinate of the micro-volume
/// \param y         logical y-coordinate of the micro-volume
/// \param faceType  type of the volume (two types of micro-faces exist)
/// \param dof       DoF ID (there may be more than one DoF per volume)
/// \param ndofs     number of DoFs per micro-volume
/// \param width     number of micro-edges per edge of the macro-face
/// \param memLayout specifies the memory layout for which the array index is computed
///
/// \return array index
inline constexpr uint_t index( uint_t                x,
                               uint_t                y,
                               facedof::FaceType     faceType,
                               uint_t                dof,
                               uint_t                ndofs,
                               uint_t                width,
                               VolumeDoFMemoryLayout memLayout )
{
   const auto numMicroVolumes = levelinfo::num_microfaces_per_face_from_width( width );

   const auto microVolume = faceType == facedof::FaceType::GRAY ? hyteg::indexing::macroFaceIndex( width, x, y ) :
                                                                  hyteg::indexing::macroFaceIndex( width - 1, x, y );

   if ( memLayout == VolumeDoFMemoryLayout::SoA )
   {
      const auto idx = numMicroVolumes * dof + microVolume;
      return idx;
   }
   else
   {
      const auto idx = microVolume * ndofs + dof;
      return idx;
   }
}

/// \brief Returns the array-index of the specified 2D volume DoF on a ghost-layer.
///
/// \param x         logical x-coordinate of the micro-volume on the ghost layer
/// \param dof       DoF ID (there may be more than one DoF per volume)
/// \param ndofs     number of DoFs per micro-volume
/// \param width     number of micro-edges per edge of the macro-face
/// \param memLayout specifies the memory layout for which the array index is computed
///
/// \return array index
inline constexpr uint_t indexGhostLayer( uint_t x, uint_t dof, uint_t ndofs, uint_t width, VolumeDoFMemoryLayout memLayout )
{
   const auto numMicroVolumes = width;

   const auto microVolume = hyteg::indexing::macroEdgeIndex( width, x );

   if ( memLayout == VolumeDoFMemoryLayout::SoA )
   {
      const auto idx = numMicroVolumes * dof + microVolume;
      return idx;
   }
   else
   {
      const auto idx = microVolume * ndofs + dof;
      return idx;
   }
}

} // namespace indexing

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
                      VolumeDoFMemoryLayout                      memoryLayout );

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
                      VolumeDoFMemoryLayout                      memoryLayout );

//   /// \brief Sets all DoFs to a specified constant.
//   void setConstant( real_t constant, uint_t level );
//
//   /// \brief Sets all DoFs to 0.
//   void setZero( uint_t level );
//
//   /// \brief Assigns a linear combination of multiple VolumeDoFFunctions to this.
//   void assign( const std::vector< ValueType >&                                                      scalars,
//                const std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > >& functions,
//                uint_t                                                                               level );

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
      return data[indexing::index( idx.x(),
                                   idx.y(),
                                   faceType,
                                   dofID,
                                   numScalarsPerPrimitive_.at( primitiveID ),
                                   indexing::levelToWidth( level ),
                                   memoryLayout_ )];
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
      return data[indexing::index( idx.x(),
                                   idx.y(),
                                   faceType,
                                   dofID,
                                   numScalarsPerPrimitive_.at( primitiveID ),
                                   indexing::levelToWidth( level ),
                                   memoryLayout_ )];
   }

 private:
   void allocateMemory();

   std::string name_;

   std::shared_ptr< PrimitiveStorage > storage_;

   uint_t minLevel_;
   uint_t maxLevel_;

   std::map< PrimitiveID, uint_t > numScalarsPerPrimitive_;

   VolumeDoFMemoryLayout memoryLayout_;

   /// Each data ID stores the inner scalar fields for all levels.
   PrimitiveDataID< FunctionMemory< ValueType >, Face > faceInnerDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell > cellInnerDataID_;

   /// One data ID per ghost-layer. Should be up to 3 in 2D and 4 in 3D.
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > faceGhostLayerDataIDs_;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > cellGhostLayerDataIDs_;
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