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

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

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
    index( idx_t x, idx_t y, facedof::FaceType faceType, uint_t dof, uint_t ndofs, uint_t level, VolumeDoFMemoryLayout memLayout )
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

/// \brief Returns the array-index of the specified 3D volume DoF.
///
/// \param x         logical x-coordinate of the micro-volume
/// \param y         logical y-coordinate of the micro-volume
/// \param z         logical z-coordinate of the micro-volume
/// \param cellType  type of the volume (six types of micro-cells exist)
/// \param dof       DoF ID (there may be more than one DoF per volume)
/// \param ndofs     number of DoFs per micro-volume
/// \param level     refinement level
/// \param memLayout specifies the memory layout for which the array index is computed
///
/// \return array index
inline constexpr uint_t index( idx_t                 x,
                               idx_t                 y,
                               idx_t                 z,
                               celldof::CellType     cellType,
                               uint_t                dof,
                               uint_t                ndofs,
                               uint_t                level,
                               VolumeDoFMemoryLayout memLayout )
{
   const auto numMicroVolumes = levelinfo::num_microcells_per_cell( level );

   const auto microVolume = celldof::macrocell::index( level, x, y, z, cellType );

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
/// This function is used to give direct access to the ghost-layers via ghost-layer local indexing.
///
/// \param x           logical x-coordinate of the micro-volume
/// \param dof         DoF ID (there may be more than one DoF per volume)
/// \param ndofs       number of DoFs per micro-volume on the respective neighbor primitive
/// \param level       refinement level
/// \param memLayout   specifies the memory layout for which the array index is computed
///
/// \return array index
inline uint_t indexGhostLayerDirectly( idx_t x, uint_t dof, uint_t ndofs, uint_t level, VolumeDoFMemoryLayout memLayout )
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
/// \brief Returns the array-index of the specified 3D volume DoF on a ghost-layer.
///
/// This function is used to give direct access to the ghost-layers via ghost-layer local indexing.
///
/// Notes on cell-types: the cell type cannot be white-down, since those elements do not share a facet with the macro-cell
/// boundary. All types that are blue-* or green-* yield the same array indices.
///
/// \param x           logical x-coordinate of the micro-volume
/// \param y           logical y-coordinate of the micro-volume
/// \param cellType    local cell type of the micro-volume
/// \param dof         DoF ID (there may be more than one DoF per volume)
/// \param ndofs       number of DoFs per micro-volume on the respective neighbor primitive
/// \param level       refinement level
/// \param memLayout   specifies the memory layout for which the array index is computed
///
/// \return array index
inline uint_t indexGhostLayerDirectly( idx_t                 x,
                                       idx_t                 y,
                                       celldof::CellType     cellType,
                                       uint_t                dof,
                                       uint_t                ndofs,
                                       uint_t                level,
                                       VolumeDoFMemoryLayout memLayout )
{
   WALBERLA_ASSERT_NOT_IDENTICAL( cellType, celldof::CellType::WHITE_DOWN, "Invalid cell type." );

   const auto numMicroVolumes = levelinfo::num_microfaces_per_face( level );

   uint_t microVolume;

   if ( cellType == celldof::CellType::WHITE_UP )
   {
      microVolume = facedof::macroface::index( level, x, y, facedof::FaceType::GRAY );
   }
   else
   {
      microVolume = facedof::macroface::index( level, x, y, facedof::FaceType::BLUE );
   }

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

/// \brief Given a logical micro-volume index of a micro-volume that has contact to the macro-volume boundary, this function
///        returns the corresponding neighboring micro-cell index in the ghost-layer data structure.
///
/// This function exploits that there is only exactly one neighboring micro-volume per ghost-layer and micro-volume.
///
/// \param localEdgeID   the local macro-interface ID of the corresponding boundary
/// \param xInner        logical x index of the (inner) micro-element whose neighbor we want to access
/// \param yInner        logical y index of the (inner) micro-element whose neighbor we want to access
/// \param faceTypeInner element type of the (inner) micro-element whose neighbor we want to access
/// \param dof           dof index
/// \param ndofs         number of dofs in the ghost layer
/// \param level         refinement level
/// \param memLayout     memory layout
/// \return array index for access in the ghost-layer
inline uint_t indexNeighborInGhostLayer( uint_t                localEdgeID,
                                         idx_t                 xInner,
                                         idx_t                 yInner,
                                         facedof::FaceType     faceTypeInner,
                                         uint_t                dof,
                                         uint_t                ndofs,
                                         uint_t                level,
                                         VolumeDoFMemoryLayout memLayout )
{
   WALBERLA_ASSERT_EQUAL( faceTypeInner, facedof::FaceType::GRAY, "All volumes at the boundary are of face type GRAY." );
   WALBERLA_UNUSED( faceTypeInner );

   switch ( localEdgeID )
   {
   case 0:
      // bottom edge
      WALBERLA_ASSERT_EQUAL( yInner, 0, "Iterating over bottom edge, yInner must be 0." );
      WALBERLA_UNUSED( yInner );
      return indexGhostLayerDirectly( xInner, dof, ndofs, level, memLayout );
   case 1:
      // left edge
      WALBERLA_ASSERT_EQUAL( xInner, 0, "Iterating over left edge, xInner must be 0." );
      WALBERLA_UNUSED( xInner );
      return indexGhostLayerDirectly( yInner, dof, ndofs, level, memLayout );
   case 2:
      // diagonal edge
      WALBERLA_ASSERT_EQUAL( xInner + yInner,
                             levelinfo::num_microedges_per_edge( level ) - 1,
                             "Iterating over diagonal edge, xInner + yInner must be "
                                 << levelinfo::num_microedges_per_edge( level ) - 1 << "." );
      return indexGhostLayerDirectly( yInner, dof, ndofs, level, memLayout );
   default:
      WALBERLA_ABORT( "Invalid local macro-edge ID: " << localEdgeID );
   }
}

/// \brief Given a logical micro-volume index of a micro-volume that has contact to the macro-volume boundary, this function
///        returns the corresponding neighboring micro-cell index in the ghost-layer data structure.
///
/// This function exploits that there is only exactly one neighboring micro-volume per ghost-layer and micro-volume.
///
/// \param localEdgeID   the local macro-interface ID of the corresponding boundary
/// \param xInner        logical x index of the (inner) micro-element whose neighbor we want to access
/// \param yInner        logical y index of the (inner) micro-element whose neighbor we want to access
/// \param zInner        logical z index of the (inner) micro-element whose neighbor we want to access
/// \param cellTypeInner element type of the (inner) micro-element whose neighbor we want to access
/// \param dof           dof index
/// \param ndofs         number of dofs in the ghost layer
/// \param level         refinement level
/// \param memLayout     memory layout
/// \return array index for access in the ghost-layer
inline uint_t indexNeighborInGhostLayer( uint_t                localFaceID,
                                         idx_t                 xInner,
                                         idx_t                 yInner,
                                         idx_t                 zInner,
                                         celldof::CellType     cellTypeInner,
                                         uint_t                dof,
                                         uint_t                ndofs,
                                         uint_t                level,
                                         VolumeDoFMemoryLayout memLayout )
{
   WALBERLA_ASSERT_UNEQUAL( cellTypeInner, celldof::CellType::WHITE_DOWN, "Cell type WHITE_DOWN has no boundary contact." );

   switch ( localFaceID )
   {
   case 0:
      // bottom face
      WALBERLA_ASSERT_EQUAL( zInner, 0, "Iterating over bottom face, zInner must be 0." );
      WALBERLA_UNUSED( zInner );
      return indexGhostLayerDirectly( xInner, yInner, cellTypeInner, dof, ndofs, level, memLayout );
   case 1:
      // front face
      WALBERLA_ASSERT_EQUAL( yInner, 0, "Iterating over front face, yInner must be 0." );
      WALBERLA_UNUSED( yInner );
      return indexGhostLayerDirectly( xInner, zInner, cellTypeInner, dof, ndofs, level, memLayout );
   case 2:
      // left face
      WALBERLA_ASSERT_EQUAL( xInner, 0, "Iterating over left face, xInner must be 0." );
      WALBERLA_UNUSED( xInner );
      return indexGhostLayerDirectly( yInner, zInner, cellTypeInner, dof, ndofs, level, memLayout );
   case 3:
      // back face
      WALBERLA_DEBUG_SECTION()
      {
         if ( cellTypeInner == celldof::CellType::WHITE_UP )
         {
            WALBERLA_ASSERT_EQUAL( xInner + yInner + zInner,
                                   levelinfo::num_microedges_per_edge( level ) - 1,
                                   "Iterating over back face, cell type "
                                       << celldof::CellTypeToStr.at( cellTypeInner ) << ", sum of all index coordinates must be "
                                       << levelinfo::num_microedges_per_edge( level ) - 1 << "." );
         }
         else
         {
            WALBERLA_ASSERT_EQUAL( xInner + yInner + zInner,
                                   levelinfo::num_microedges_per_edge( level ) - 2,
                                   "Iterating over back face, cell type "
                                       << celldof::CellTypeToStr.at( cellTypeInner ) << ", sum of all index coordinates must be "
                                       << levelinfo::num_microedges_per_edge( level ) - 1 << "." );
         }
      }
      WALBERLA_UNUSED( xInner );
      return indexGhostLayerDirectly( yInner, zInner, cellTypeInner, dof, ndofs, level, memLayout );
   default:
      WALBERLA_ABORT( "Invalid local macro-face ID: " << localFaceID );
   }
}

///@}

/// \brief Given a micro-element on a refined macro element, this function returns the "parent"-micro-element on the next coarser
/// macro.
inline void getCoarseMicroElementFromFineMicroElement( const hyteg::indexing::Index& fineElementIdx,
                                                       const facedof::FaceType&      fineFaceType,
                                                       hyteg::indexing::Index&       coarseElementIdx,
                                                       facedof::FaceType&            coarseFaceType )
{
   const auto x = fineElementIdx.x();
   const auto y = fineElementIdx.y();

   if ( fineFaceType == facedof::FaceType::GRAY )
   {
      coarseElementIdx.x() = x / 2;
      coarseElementIdx.y() = y / 2;

      if ( x % 2 != 0 && y % 2 != 0 )
      {
         coarseFaceType = facedof::FaceType::BLUE;
      }
      else
      {
         coarseFaceType = fineFaceType;
      }
   }
   else
   {
      if ( x % 2 == 0 && y % 2 == 0 )
      {
         coarseFaceType = facedof::FaceType::GRAY;
      }
      else
      {
         coarseFaceType = fineFaceType;
      }
   }
}

/// \brief Given a micro-element on a coarse macro element, this function returns the "child"-micro-elements on the next finer
/// macro.
inline void getFineMicroElementsFromCoarseMicroElement( const hyteg::indexing::Index&          coarseElementIdx,
                                                        const facedof::FaceType&               coarseFaceType,
                                                        std::vector< hyteg::indexing::Index >& fineElementIdx,
                                                        std::vector< facedof::FaceType >&      fineFaceType )
{
   fineElementIdx.clear();
   fineFaceType.clear();

   fineElementIdx.resize( 4 );
   fineFaceType.resize( 4 );

   if ( coarseFaceType == facedof::FaceType::GRAY )
   {
      fineElementIdx[0] = 2 * coarseElementIdx;
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineElementIdx[3] = 2 * coarseElementIdx;

      fineFaceType[0] = facedof::FaceType::GRAY;
      fineFaceType[1] = facedof::FaceType::GRAY;
      fineFaceType[2] = facedof::FaceType::GRAY;
      fineFaceType[3] = facedof::FaceType::BLUE;
   }
   else
   {
      fineElementIdx[0] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineElementIdx[3] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );

      fineFaceType[0] = facedof::FaceType::BLUE;
      fineFaceType[1] = facedof::FaceType::BLUE;
      fineFaceType[2] = facedof::FaceType::BLUE;
      fineFaceType[3] = facedof::FaceType::GRAY;
   }
}

/// \brief Given a micro-element on a coarse macro element, this function returns the "child"-micro-elements on the next finer
/// macro.
inline void getFineMicroElementsFromCoarseMicroElement( const hyteg::indexing::Index&          coarseElementIdx,
                                                        const celldof::CellType&               coarseCellType,
                                                        std::vector< hyteg::indexing::Index >& fineElementIdx,
                                                        std::vector< celldof::CellType >&      fineCellType )
{
   fineElementIdx.clear();
   fineCellType.clear();

   fineElementIdx.resize( 8 );
   fineCellType.resize( 8 );

   if ( coarseCellType == celldof::CellType::WHITE_UP )
   {
      fineElementIdx[0] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 0 );
      fineCellType[0]   = celldof::CellType::WHITE_UP;
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineCellType[1]   = celldof::CellType::WHITE_UP;
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineCellType[2]   = celldof::CellType::WHITE_UP;
      fineElementIdx[3] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 1 );
      fineCellType[3]   = celldof::CellType::WHITE_UP;
      fineElementIdx[4] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 0 );
      fineCellType[4]   = celldof::CellType::BLUE_UP;
      fineElementIdx[5] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 0 );
      fineCellType[5]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[6] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 0 );
      fineCellType[6]   = celldof::CellType::GREEN_UP;
      fineElementIdx[7] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 0 );
      fineCellType[7]   = celldof::CellType::GREEN_DOWN;
   }
   else if ( coarseCellType == celldof::CellType::WHITE_DOWN )
   {
      fineElementIdx[0] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );
      fineCellType[0]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 1 );
      fineCellType[1]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 1 );
      fineCellType[2]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[3] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 1 );
      fineCellType[3]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[4] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 1 );
      fineCellType[4]   = celldof::CellType::BLUE_UP;
      fineElementIdx[5] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 1 );
      fineCellType[5]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[6] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 1 );
      fineCellType[6]   = celldof::CellType::GREEN_UP;
      fineElementIdx[7] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 1 );
      fineCellType[7]   = celldof::CellType::GREEN_DOWN;
   }
   else if ( coarseCellType == celldof::CellType::BLUE_UP )
   {
      fineElementIdx[0] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );
      fineCellType[0]   = celldof::CellType::WHITE_UP;
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineCellType[1]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineCellType[2]   = celldof::CellType::BLUE_UP;
      fineElementIdx[3] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineCellType[3]   = celldof::CellType::BLUE_UP;
      fineElementIdx[4] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );
      fineCellType[4]   = celldof::CellType::BLUE_UP;
      fineElementIdx[5] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 1 );
      fineCellType[5]   = celldof::CellType::BLUE_UP;
      fineElementIdx[6] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );
      fineCellType[6]   = celldof::CellType::GREEN_UP;
      fineElementIdx[7] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineCellType[7]   = celldof::CellType::GREEN_DOWN;
   }
   else if ( coarseCellType == celldof::CellType::BLUE_DOWN )
   {
      fineElementIdx[0] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 1 );
      fineCellType[0]   = celldof::CellType::WHITE_UP;
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 1 );
      fineCellType[1]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineCellType[2]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[3] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 1 );
      fineCellType[3]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[4] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 1 );
      fineCellType[4]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[5] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 1 );
      fineCellType[5]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[6] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 1 );
      fineCellType[6]   = celldof::CellType::GREEN_UP;
      fineElementIdx[7] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 1 );
      fineCellType[7]   = celldof::CellType::GREEN_DOWN;
   }
   else if ( coarseCellType == celldof::CellType::GREEN_UP )
   {
      fineElementIdx[0] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 1 );
      fineCellType[0]   = celldof::CellType::WHITE_UP;
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 0 );
      fineCellType[1]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 1 );
      fineCellType[2]   = celldof::CellType::BLUE_UP;
      fineElementIdx[3] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineCellType[3]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[4] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 0 );
      fineCellType[4]   = celldof::CellType::GREEN_UP;
      fineElementIdx[5] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineCellType[5]   = celldof::CellType::GREEN_UP;
      fineElementIdx[6] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 0, 1 );
      fineCellType[6]   = celldof::CellType::GREEN_UP;
      fineElementIdx[7] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 1 );
      fineCellType[7]   = celldof::CellType::GREEN_UP;
   }
   else
   {
      fineElementIdx[0] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 1 );
      fineCellType[0]   = celldof::CellType::WHITE_UP;
      fineElementIdx[1] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineCellType[1]   = celldof::CellType::WHITE_DOWN;
      fineElementIdx[2] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 1 );
      fineCellType[2]   = celldof::CellType::BLUE_UP;
      fineElementIdx[3] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );
      fineCellType[3]   = celldof::CellType::BLUE_DOWN;
      fineElementIdx[4] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 0 );
      fineCellType[4]   = celldof::CellType::GREEN_DOWN;
      fineElementIdx[5] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 1, 0 );
      fineCellType[5]   = celldof::CellType::GREEN_DOWN;
      fineElementIdx[6] = 2 * coarseElementIdx + hyteg::indexing::Index( 1, 0, 1 );
      fineCellType[6]   = celldof::CellType::GREEN_DOWN;
      fineElementIdx[7] = 2 * coarseElementIdx + hyteg::indexing::Index( 0, 1, 1 );
      fineCellType[7]   = celldof::CellType::GREEN_DOWN;
   }
}

/// \brief Simple class that computes and holds information about the neighborhood of a volume element.
class ElementNeighborInfo
{
 public:
   using Point    = Eigen::Matrix< real_t, 3, 1 >;
   using Index    = hyteg::indexing::Index;
   using FaceType = facedof::FaceType;
   using CellType = celldof::CellType;

   ElementNeighborInfo() = default;

   ElementNeighborInfo( Index                                      elementIdx,
                        FaceType                                   faceType,
                        uint_t                                     level,
                        BoundaryCondition                          boundaryCondition,
                        PrimitiveID                                faceID,
                        const std::shared_ptr< PrimitiveStorage >& storage );

   ElementNeighborInfo( Index                                      elementIdx,
                        CellType                                   cellType,
                        uint_t                                     level,
                        BoundaryCondition                          boundaryCondition,
                        PrimitiveID                                cellID,
                        const std::shared_ptr< PrimitiveStorage >& storage );

   const Index& elementIdx() const { return elementIdx_; }

   const FaceType& faceType() const { return faceType_; }
   const CellType& cellType() const { return cellType_; }

   [[nodiscard]] const std::vector< Point >& elementVertexCoords() const { return vertexCoordsVolume_; }

   [[nodiscard]] const std::vector< Index >& elementVertexIndices() const { return vertexIndicesVolume_; }

   Index neighborElementIndices( uint_t neighbor ) const { return neighborElementIndices_[neighbor]; }

   [[nodiscard]] const std::vector< Point >& neighborElementVertexCoords( uint_t neighbor ) const
   {
      return neighborElementVertexCoords_[neighbor];
   }

   /// \brief Returns a new ElementNeighborInfo instance with updated neighborhood information for the case that the neighboring
   ///        micro-element is located on a different (neighboring) macro-volume.
   ElementNeighborInfo updateForMacroBoundary( uint_t neighbor ) const;

   [[nodiscard]] const std::vector< Point >& interfaceVertexCoords( uint_t neighbor ) const
   {
      return interfaceVertexCoords_[neighbor];
   }

   [[nodiscard]] const std::vector< Index >& interfaceVertexIndices( uint_t neighbor ) const
   {
      return interfaceVertexIndices_[neighbor];
   }

   [[nodiscard]] const Point& oppositeVertexCoords( uint_t neighbor ) const { return oppositeVertexCoords_[neighbor]; }

   [[nodiscard]] const Point& neighborOppositeVertexCoords( uint_t neighbor ) const
   {
      return neighborOppositeVertexCoords_[neighbor];
   }

   [[nodiscard]] const Point& outwardNormal( uint_t neighbor ) const { return outwardNormal_[neighbor]; }

   FaceType neighborFaceType( uint_t neighbor ) const { return neighborFaceElementTypes_[neighbor]; }

   CellType neighborCellType( uint_t neighbor ) const { return neighborCellElementTypes_[neighbor]; }

   bool   atMacroBoundary( uint_t neighbor ) const { return macroBoundary_[neighbor] != NOT_AT_BOUNDARY; }
   uint_t macroBoundaryID( uint_t neighbor ) const { return macroBoundary_[neighbor]; }

   DoFType neighborBoundaryType( uint_t neighbor ) const { return neighborBoundaryType_[neighbor]; }

 private:
   /// Dimensionality of the volume element.
   int dim_;

   /// Logical index of the element.
   Index elementIdx_;

   /// Primitive ID of the containing volume primitive.
   PrimitiveID volumeID_;

   std::shared_ptr< PrimitiveStorage > storage_;

   /// Type of the element (if 2D)
   FaceType faceType_;

   /// Type of the element (if 3D)
   CellType cellType_;

   /// Refinement level.
   uint_t level_;

   /// Logical indices of the elements' vertices.
   std::vector< Index > vertexIndicesVolume_;

   /// Coordinates of the elements' vertices.
   std::vector< Point > vertexCoordsVolume_;

   /// Logical indices of the neighbor elements.
   std::vector< Index > neighborElementIndices_;

   /// Types of the neighbor elements (if 2D).
   std::vector< FaceType > neighborFaceElementTypes_;

   /// Types of the neighbor elements (if 3D).
   std::vector< CellType > neighborCellElementTypes_;

   /// Coordinates of the neighboring elements' vertices.
   std::vector< std::vector< Point > > neighborElementVertexCoords_;

   /// Logical vertex indices of the interfaces to the neighboring elements.
   std::vector< std::vector< Index > > interfaceVertexIndices_;

   /// Coordinates of the vertices at the interfaces to the neighboring elements.
   std::vector< std::vector< Point > > interfaceVertexCoords_;

   /// Normal to the interface pointing away from the element.
   std::vector< Point > outwardNormal_;

   /// Logical vertex index of the element's vertex that is not on the interface.
   std::vector< Index > oppositeVertexIndex_;

   /// Coordinates of the element's vertex that is not on the interface.
   std::vector< Point > oppositeVertexCoords_;

   /// Logical vertex index of the neighboring element's vertex that is not on the interface.
   std::vector< Index > neighborOppositeVertexIndex_;

   /// Coordinates of the neighboring element's vertex that is not on the interface.
   std::vector< Point > neighborOppositeVertexCoords_;

   const static uint_t NOT_AT_BOUNDARY;
   /// Stores for each micro-element interface which macro-local boundary it touches.
   /// If it is not located at a boundary, the value is set to NOT_AT_BOUNDARY.
   std::vector< uint_t > macroBoundary_;

   /// Stores the DoFType (think boundary condition) for each element interface.
   std::vector< DoFType > neighborBoundaryType_;
};

// TODO: move some of these functions to some other file(s) ...

/// \brief Given a micro-vertex idx on a macro (micro ref. level l+1), this function returns the corresponding micro-vertex idx
///        on a refined macro (micro ref. level l).
///
/// Since a micro-vertex may lie on multiple finer macros, the corresponding finer macro must be specified.
/// The indices of the local macro-vertices are assumed to be chosen as in Bey (1995): Tetrahedral Grid Refinement.
///
/// Note that this function also works in 2D - all z-coordinates must simply be set to zero. Note that GRAY faces correspond
/// to WHITE_UP cells and BLUE faces to BLUE_UP cells.
///
/// This function does not check whether the returned index is inside the fine macro.
///
/// \param idxCoarse     micro-vertex idx on the coarse tet with refinement level = levelFine + 1
/// \param levelFine     refinement level of the fine macro
/// \param macroCellType cell type of the fine macro (a macro cell is refined into eight tets, four of which are of type WHITE_UP,
///                      the others are BLUE_UP, BLUE_DOWN, GREEN_UP, GREEN_DOWN
/// \param macroCellIdx  the macro cell of choice if cell type is WHITE_UP - since there are four (numbered in xyz-ordering)
/// \return the index local to the refined macro cell
inline hyteg::indexing::Index getMicroVertexIdxOnRefinedMacro( const hyteg::indexing::Index& idxCoarse,
                                                               const uint_t&                 levelFine,
                                                               const celldof::CellType&      macroCellType,
                                                               const uint_t&                 macroCellIdx )
{
   // The given level is the refinement level of the fine macros.
   // We assume that the refinement level on the coarse macro is level + 1 to match the grid structure.
   //
   // To find the corresponding micro idx on the given refined macro we exploit the fixed ordering of the local vertices of the
   // refined macros according to Bey.
   //
   // The procedure is as follows: first we determine the micro-vertex idx b on the coarse macro that corresponds to the local
   // (0, 0, 0) micro-vertex idx on the respective fine macro.
   //
   // We calculate the distance d = x - b (where x is the given idx on the coarse macro) that corresponds to the offset
   // from the fine-macro local (0, 0, 0) idx in the coordinate system of the coarse macro.
   //
   // It remains to translate this to the fine-macro coordinate system. Let f1, f2, f3 in Z^3 be the three x, y, z idx increments
   // in the coarse-macro coordinate system. Then we want to find the fine-macro local idx y in Z^3 so that:
   //
   //     y1 * f1 + y2 * f2 + y3 * f3 = d
   //
   // Eventually we need to invert the 3x3 system (f1 f2 f3). We can precompute these systems analytically, and thus only need to
   // apply the inverse. For the WHITE_UP cells with the Bey vertex-ordering, these matrices are just identities btw.

   using hyteg::celldof::CellType;
   using hyteg::indexing::Index;

   const auto widthFine = idx_t( levelinfo::num_microvertices_per_edge( levelFine ) );

   // Base idx, the three rows of the inverse of (f1 f2 f3), and the result.
   Index b, fRow1, fRow2, fRow3, idxFine;

   switch ( macroCellType )
   {
   case CellType::WHITE_UP:
      switch ( macroCellIdx )
      {
      case 0:
         b = Index( 0, 0, 0 );
         break;
      case 1:
         b = Index( widthFine - 1, 0, 0 );
         break;
      case 2:
         b = Index( 0, widthFine - 1, 0 );
         break;
      case 3:
         b = Index( 0, 0, widthFine - 1 );
         break;
      default:
         WALBERLA_ABORT( "Invalid macro cell idx." );
      }
      fRow1 = Index( 1, 0, 0 );
      fRow2 = Index( 0, 1, 0 );
      fRow3 = Index( 0, 0, 1 );
      break;
   case CellType::BLUE_UP:
      b = Index( widthFine - 1, 0, 0 );
      // (f1 f2 f3) = ( -1  0  0 )
      //              (  1  1  0 )
      //              (  0  0  1 )
      fRow1 = Index( -1, 0, 0 );
      fRow2 = Index( 1, 1, 0 );
      fRow3 = Index( 0, 0, 1 );
      // (yep the matrix is involutory :) )
      break;
   case CellType::BLUE_DOWN:
      b = Index( 0, widthFine - 1, 0 );
      // (f1 f2 f3) = (  0  1  0 )
      //              ( -1 -1  0 )
      //              (  1  1  1 )
      fRow1 = Index( -1, -1, 0 );
      fRow2 = Index( 1, 0, 0 );
      fRow3 = Index( 0, 1, 1 );
      break;
   case CellType::GREEN_UP:
      b = Index( widthFine - 1, 0, 0 );
      // (f1 f2 f3) = ( -1 -1  0 )
      //              (  1  0  0 )
      //              (  0  1  1 )
      fRow1 = Index( 0, 1, 0 );
      fRow2 = Index( -1, -1, 0 );
      fRow3 = Index( 1, 1, 1 );
      break;
   case CellType::GREEN_DOWN:
      b = Index( 0, widthFine - 1, 0 );
      // (f1 f2 f3) = (  1  1  0 )
      //              (  0 -1  0 )
      //              (  0  1  1 )
      fRow1 = Index( 1, 1, 0 );
      fRow2 = Index( 0, -1, 0 );
      fRow3 = Index( 0, 1, 1 );
      // (yep the matrix is involutory, too :) )
      break;
   }

   auto d = idxCoarse - b;

   idxFine[0] = fRow1.dot( d );
   idxFine[1] = fRow2.dot( d );
   idxFine[2] = fRow3.dot( d );

   return idxFine;
}

/// \brief Given a micro-vertex idx on a macro (micro ref. level l), this function returns the corresponding micro-vertex idx
///        on the next coarser macro (micro ref. level l + 1).
///
/// The indices of the local macro-vertices are assumed to be chosen as in Bey (1995): Tetrahedral Grid Refinement.
///
/// Note that this function also works in 2D - all z-coordinates must simply be set to zero. Note that GRAY faces correspond
/// to WHITE_UP cells and BLUE faces to BLUE_UP cells.
///
/// This function does not check whether the returned index is inside the coarse macro.
///
/// \param idxFine       micro-vertex idx on the fine tet with refinement level = levelFine
/// \param levelFine     refinement level of the fine macro
/// \param macroCellType cell type of the fine macro (a macro cell is refined into eight tets, four of which are of type WHITE_UP,
///                      the others are BLUE_UP, BLUE_DOWN, GREEN_UP, GREEN_DOWN
/// \param macroCellIdx  the macro cell of choice if cell type is WHITE_UP - since there are four (numbered in xyz-ordering)
/// \return the index local to the refined macro cell
inline hyteg::indexing::Index getMicroVertexIdxOnCoarserMacro( const hyteg::indexing::Index& idxFine,
                                                               const uint_t&                 levelFine,
                                                               const celldof::CellType&      macroCellType,
                                                               const uint_t&                 macroCellIdx )
{
   // The given level is the refinement level of the fine macros.
   // We assume that the refinement level on the coarse macro is level + 1 to match the grid structure.
   //
   // To find the corresponding micro idx on the coarser macro we exploit the fixed ordering of the local vertices of the
   // refined macros according to Bey.
   //
   // The procedure is as follows: we first translate to coarse-macro local coordinates. But still the offset to the idx
   // that has coords (0, 0, 0) in the fine-macro space remains. So we need to add that offset to the translated index.
   //
   // This basically reverses getMicroVertexIdxOnRefinedMacro() so have a look at that implementation, too.

   using hyteg::celldof::CellType;
   using hyteg::indexing::Index;

   const auto widthFine = idx_t( levelinfo::num_microvertices_per_edge( levelFine ) );

   // Base idx, the three rows of (f1 f2 f3), and the result.
   Index b, fRow1, fRow2, fRow3, idxCoarse;

   switch ( macroCellType )
   {
   case CellType::WHITE_UP:
      switch ( macroCellIdx )
      {
      case 0:
         b = Index( 0, 0, 0 );
         break;
      case 1:
         b = Index( widthFine - 1, 0, 0 );
         break;
      case 2:
         b = Index( 0, widthFine - 1, 0 );
         break;
      case 3:
         b = Index( 0, 0, widthFine - 1 );
         break;
      default:
         WALBERLA_ABORT( "Invalid macro cell idx." );
      }
      fRow1 = Index( 1, 0, 0 );
      fRow2 = Index( 0, 1, 0 );
      fRow3 = Index( 0, 0, 1 );
      break;
   case CellType::BLUE_UP:
      b     = Index( widthFine - 1, 0, 0 );
      fRow1 = Index( -1, 0, 0 );
      fRow2 = Index( 1, 1, 0 );
      fRow3 = Index( 0, 0, 1 );
      break;
   case CellType::BLUE_DOWN:
      b     = Index( 0, widthFine - 1, 0 );
      fRow1 = Index( 0, 1, 0 );
      fRow2 = Index( -1, -1, 0 );
      fRow3 = Index( 1, 1, 1 );
      break;
   case CellType::GREEN_UP:
      b     = Index( widthFine - 1, 0, 0 );
      fRow1 = Index( -1, -1, 0 );
      fRow2 = Index( 1, 0, 0 );
      fRow3 = Index( 0, 1, 1 );
      break;
   case CellType::GREEN_DOWN:
      b     = Index( 0, widthFine - 1, 0 );
      fRow1 = Index( 1, 1, 0 );
      fRow2 = Index( 0, -1, 0 );
      fRow3 = Index( 0, 1, 1 );
      break;
   }

   idxCoarse[0] = fRow1.dot( idxFine );
   idxCoarse[1] = fRow2.dot( idxFine );
   idxCoarse[2] = fRow3.dot( idxFine );

   idxCoarse += b;

   return idxCoarse;
}

/// \brief Given a micro-volume index on level l+1, this function returns the corresponding local micro-index on the
/// corresponding refined macro.
///
/// The prerequisite is that the passed macro-volume has already been refined once.
///
/// It is relevant, that the micro-volumes of a refined macro-volume have the same shape regardless of whether
///  - the macro-volume was refined or
///  - the micro-refinement level was increased.
///
/// Now, this function receives a refinement level, and a micro-volume idx on level l+1.
/// It first computes the (refined) macro-volume that contains that index.
/// Then it computes the micro-volume idx _local_ to the refined macro-volume.
///
/// \param face            [in] coarse macro-primitive
/// \param level           [in] refinement level of the refined macro-primitives
/// \param idx             [in] micro-volume index on the coarse macro-primitive with (level + 1)
/// \param faceType        [in] micro-volume type
/// \param refinedMacroPID [out] PrimitiveID of the fine macro
/// \param outIdx          [out] micro-volume index (on the passed level) local to the fine macro
/// \param outFaceType     [out] micro-volume type local to the fine macro
inline void getVolumeIdxOnRefinedMacro( const Face&                   face,
                                        uint_t                        level,
                                        const hyteg::indexing::Index& idx,
                                        const facedof::FaceType&      faceType,
                                        PrimitiveID&                  refinedMacroPID,
                                        hyteg::indexing::Index&       outIdx,
                                        facedof::FaceType&            outFaceType )
{
   // Two steps: first we need to determine the fine macro that we end up in.
   // Then we need to translate to fine-macro local coordinates.
   //
   // The first step is performed "manually". For the second we simply translate all the micro-vertex coords of the cell and
   // then retrieve the corresponding micro-volume from those. This is a bit easier. The functions for that exist already.

   using hyteg::indexing::Index;

   WALBERLA_ASSERT( face.hasChildren(), "The passed macro-face must have children. Otherwise this function has no job." );

   uint_t widthFineMacro = levelinfo::num_microedges_per_edge( level );

   // a) Finding the correct macro.

   celldof::CellType macroType = celldof::CellType::WHITE_UP;
   uint_t            macroIdx  = 42;

   if ( faceType == facedof::FaceType::GRAY )
   {
      if ( idx.x() + idx.y() < widthFineMacro )
      {
         macroIdx = 0;
      }
      else if ( idx.x() >= widthFineMacro )
      {
         macroIdx = 1;
      }
      else if ( idx.y() >= widthFineMacro )
      {
         macroIdx = 2;
      }
      else
      {
         macroType = celldof::CellType::BLUE_UP;
      }
   }
   else
   {
      if ( idx.x() + idx.y() < widthFineMacro - 1 )
      {
         macroIdx = 0;
      }
      else if ( idx.x() >= widthFineMacro )
      {
         macroIdx = 1;
      }
      else if ( idx.y() >= widthFineMacro )
      {
         macroIdx = 2;
      }
      else
      {
         macroType = celldof::CellType::BLUE_UP;
      }
   }

   if ( macroType == celldof::CellType::WHITE_UP )
   {
      refinedMacroPID = face.childFaces().at( macroIdx );
   }
   else
   {
      refinedMacroPID = face.childFaces().at( 3 );
   }

   // b) Finding the correct micro.

   auto                   microVerticesCoarse = facedof::macroface::getMicroVerticesFromMicroFace( idx, faceType );
   std::array< Index, 3 > microVerticesFine;
   for ( uint_t i = 0; i < 3; i++ )
   {
      microVerticesFine[i] = getMicroVertexIdxOnRefinedMacro( microVerticesCoarse[i], level, macroType, macroIdx );
   }
   auto microFace = facedof::macroface::getMicroFaceFromMicroVertices( microVerticesFine );

   outIdx      = microFace.first;
   outFaceType = microFace.second;
}

/// \brief Given a micro-volume index on level l, this function returns the corresponding local micro-index on the
/// corresponding coarser macro on level l+1.
///
/// The prerequisite is that the passed macro-volume has a parent volume.
///
/// It is relevant, that the micro-volumes of a refined macro-volume have the same shape regardless of whether
///  - the macro-volume was refined or
///  - the micro-refinement level was increased.
///
/// Now, this function receives a refinement level l, and a micro-volume idx on level l.
/// It first computes the coarse macro-volume that contains that index (there is only one candidate).
/// Then it computes the micro-volume idx _local_ to the coarse macro-volume on level l + 1.
///
/// \param face            [in] fine macro-primitive
/// \param level           [in] refinement level of the refined macro-primitives
/// \param idx             [in] micro-volume index on the fine macro-primitive with level level
/// \param faceType        [in] micro-volume type
/// \param coarseMacroPID  [out] PrimitiveID of the coarse macro
/// \param outIdx          [out] micro-volume index (on level + 1) local to the coarse macro
/// \param outFaceType     [out] micro-volume type local to the coarse macro
inline void getVolumeIdxOnCoarseMacro( const PrimitiveStorage&       storage,
                                       const Face&                   face,
                                       uint_t                        level,
                                       const hyteg::indexing::Index& idx,
                                       const facedof::FaceType&      faceType,
                                       PrimitiveID&                  coarseMacroPID,
                                       hyteg::indexing::Index&       outIdx,
                                       facedof::FaceType&            outFaceType )
{
   // There is only one possible coarse macro. But we have to find out which type of macro we are on as seen from the coarse
   // macro. Then we have to translate the micro coordinates.

   using hyteg::indexing::Index;

   WALBERLA_ASSERT( face.hasParent(), "The passed macro-face must have a parent. Otherwise this function has no job." );

   coarseMacroPID = face.parent();

   const auto& coarseFace = *storage.getFace( coarseMacroPID );

   uint_t            fineMacroId;
   celldof::CellType macroType = celldof::CellType::WHITE_UP;
   for ( uint_t i = 0; i < 4; i++ )
   {
      if ( coarseFace.childFaces().at( i ) == face.getID() )
      {
         fineMacroId = i;
         break;
      }
   }

   if ( fineMacroId == 3 )
   {
      macroType = celldof::CellType::BLUE_UP;
   }

   auto                   microVerticesFine = facedof::macroface::getMicroVerticesFromMicroFace( idx, faceType );
   std::array< Index, 3 > microVerticesCoarse;
   for ( uint_t i = 0; i < 3; i++ )
   {
      microVerticesCoarse[i] = getMicroVertexIdxOnCoarserMacro( microVerticesFine[i], level, macroType, fineMacroId );
   }
   auto microFace = facedof::macroface::getMicroFaceFromMicroVertices( microVerticesCoarse );

   outIdx      = microFace.first;
   outFaceType = microFace.second;
}

} // namespace indexing
} // namespace volumedofspace
} // namespace hyteg
