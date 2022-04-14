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
#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"

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
/// \param faceTypeInner element type of the (inner) micro-element whose neighbor we want to access
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
   WALBERLA_ASSERT_UNEQUAL( cellTypeInner, celldof::CellType::WHITE_DOWN, "Cell type WHITE_DOWN has not boundary contact." );

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
                        CellType                                   faceType,
                        uint_t                                     level,
                        BoundaryCondition                          boundaryCondition,
                        PrimitiveID                                cellID,
                        const std::shared_ptr< PrimitiveStorage >& storage );

   [[nodiscard]] const std::vector< Point >& elementVertexCoords() const { return vertexCoordsVolume_; }

   Index neighborElementIndices( uint_t neighbor ) const { return neighborElementIndices_[neighbor]; }

   [[nodiscard]] const std::vector< Point >& neighborElementVertexCoords( uint_t neighbor ) const
   {
      return neighborElementVertexCoords_[neighbor];
   }

   void macroBoundaryNeighborElementVertexCoords( uint_t                neighbor,
                                                  std::vector< Point >& neighborElementVertexCoords,
                                                  Point&                neighborOppositeVertexCoords ) const;

   [[nodiscard]] const std::vector< Point >& interfaceVertexCoords( uint_t neighbor ) const
   {
      return interfaceVertexCoords_[neighbor];
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

} // namespace indexing
} // namespace volumedofspace
} // namespace hyteg