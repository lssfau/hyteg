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

namespace hyteg {
namespace volumedofspace {
namespace indexing {

/// \brief Defines the memory layout for a VolumeDoFFunction.
///
/// structure-of-arrays (SoA): the innermost variable determines the index of the current micro-volume
/// array-of-structures (AoS): the innermost variable determines the DoF of the current micro-volume
enum class VolumeDoFMemoryLayout
{
   SoA,
   AoS
};

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
/// \param level     refinement level
/// \param memLayout specifies the memory layout for which the array index is computed
///
/// \return array index
inline constexpr uint_t
    index( int x, int y, facedof::FaceType faceType, uint_t dof, uint_t ndofs, uint_t level, VolumeDoFMemoryLayout memLayout )
{
   const auto numMicroVolumes = levelinfo::num_microfaces_per_face( level );

   const auto microVolume = facedof::macroface::index( level, x, y, faceType );

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

///@{
/// @name Ghost-layer indexing functions.
///
/// There are distinct arrays for each ghost-layer of a macro-volume (3 arrays for macro-faces, 4 arrays for macro-cells).
/// Each array is indexed starting from the lowest local vertex ID to the largest.
/// In 3D, the inner-most loop-index is the direction from smallest to second largest local vertex ID.
/// Also, in 3D there are two micro-cell types per ghost-layer.
///

/// \brief Returns the array-index of the specified 2D volume DoF on a ghost-layer.
///
/// This function is used to give direct access to the ghost-layers via ghost-layer loca indexing.
///´
/// \param x           logical x-coordinate of the micro-volume
/// \param dof         DoF ID (there may be more than one DoF per volume)
/// \param ndofs       number of DoFs per micro-volume on the respective neighbor primitive
/// \param level       refinement level
/// \param memLayout   specifies the memory layout for which the array index is computed
///
/// \return array index
inline uint_t indexGhostLayerDirectly( int x, uint_t dof, uint_t ndofs, uint_t level, VolumeDoFMemoryLayout memLayout )
{
   const auto numMicroVolumes = levelinfo::num_microedges_per_edge( level );
   const auto microVolume     = hyteg::indexing::macroEdgeIndex( numMicroVolumes, x );

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

/// \brief Returns the array-index of the specified 2D volume DoF on a ghost-layer via extended indices.
///
/// To seamlessly access neighboring volumes that are located on ghost-layers indexing functions are supplied that accept negative
/// logical indices and element types. The functions only require information about which macro-interface is considered.
///
/// For example, when iterating over a 2D macro-face at the boundary macro-edge with local ID 1, accessing the left neighbor
/// of a micro-volume, the index can be calculated as follows:
///
/// \code{.cpp}
///   int microFaceX, microFaceY;
///   // fill ...
///   FaceType faceType         = FaceType::GRAY;
///   FaceType neighborFaceType = FaceType::BLUE; // always the same in 2D
///   uint_t localEdgeID        = 1;
///
///   uint_t arrIdxBottomNeighbor = indexGhostLayer( localEdgeID, microFaceX - 1, microFaceY, neighborFaceType, dof, nDofsNeighbor, level, memLayout );
/// \endcode
///
/// \param localEdgeID local ID of the neighboring macro-edge
/// \param x           logical x-coordinate of the micro-volume
/// \param y           logical y-coordinate of the micro-volume
/// \param faceType    type of the volume from a local perspective (as if the macro-primitive was larger)
/// \param dof         DoF ID (there may be more than one DoF per volume)
/// \param ndofs       number of DoFs per micro-volume on the respective neighbor primitive
/// \param level       refinement level
/// \param memLayout   specifies the memory layout for which the array index is computed
///
/// \return array index
inline uint_t indexGhostLayer( uint_t                localEdgeID,
                               int                   x,
                               int                   y,
                               facedof::FaceType     faceType,
                               uint_t                dof,
                               uint_t                ndofs,
                               uint_t                level,
                               VolumeDoFMemoryLayout memLayout )
{
   WALBERLA_ASSERT_EQUAL(
       faceType, facedof::FaceType::BLUE, "All boundary edges are of face type BLUE (from a macro-local perspective)." );

   const auto numMicroVolumes = levelinfo::num_microedges_per_edge( level );
   uint_t     microVolume     = std::numeric_limits< uint_t >::max();

   switch ( localEdgeID )
   {
   case 0:
      // bottom edge
      WALBERLA_ASSERT_EQUAL( y, -1, "Iterating over bottom edge, y must be -1." );
      WALBERLA_UNUSED( y );
      return indexGhostLayerDirectly( x, dof, ndofs, level, memLayout );
   case 1:
      // left edge
      WALBERLA_ASSERT_EQUAL( x, -1, "Iterating over left edge, x must be -1." );
      WALBERLA_UNUSED( x );
      return indexGhostLayerDirectly( y, dof, ndofs, level, memLayout );
   case 2:
      // diagonal edge
      WALBERLA_ASSERT_EQUAL(
          x + y, numMicroVolumes - 1, "Iterating over diagonal edge, x+y must be " << numMicroVolumes - 1 << "." );
      return indexGhostLayerDirectly( y, dof, ndofs, level, memLayout );
   default:
      WALBERLA_ABORT( "Invalid local macro-edge ID: " << localEdgeID );
   }
}

///@}

} // namespace indexing
} // namespace volumedofspace
} // namespace hyteg