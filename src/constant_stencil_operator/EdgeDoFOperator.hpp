/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2Elements3D.hpp"

namespace hyteg {

template < class EdgeDoFForm >
class EdgeDoFOperator final : public Operator< hyteg::EdgeDoFFunction< real_t >, hyteg::EdgeDoFFunction< real_t > >
{
 public:
   EdgeDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );
   EdgeDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                    size_t                                     minLevel,
                    size_t                                     maxLevel,
                    const EdgeDoFForm&                         form );
   ~EdgeDoFOperator() override = default;

   void apply( const EdgeDoFFunction< real_t >& src,
               const EdgeDoFFunction< real_t >& dst,
               uint_t                           level,
               DoFType                          flag,
               UpdateType                       updateType ) const override;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EdgeDoFFunction< idx_t >&             src,
                  const EdgeDoFFunction< idx_t >&             dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override;

   const PrimitiveDataID< StencilMemory< real_t >, Edge >&                             getEdgeStencilID() const;
   const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >& getEdgeStencil3DID() const;
   const PrimitiveDataID< StencilMemory< real_t >, Face >&                             getFaceStencilID() const;
   const PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >& getFaceStencil3DID() const;
   const PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >& getCellStencilID() const;

 private:
   void assembleStencils();

   PrimitiveDataID< StencilMemory< real_t >, Edge >                             edgeStencilID_;
   PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge > edgeStencil3DID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >                             faceStencilID_;
   PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face > faceStencil3DID_;
   PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell > cellStencilID_;

   EdgeDoFForm form_;
};

/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive );

/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive );

uint_t macroCellEdgeDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive );

template < typename EdgeDoFForm >
void assembleEdgeToEdgeStencils(
    const std::shared_ptr< PrimitiveStorage >&                                          storage,
    const uint_t&                                                                       minLevel,
    const uint_t&                                                                       maxLevel,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >& macroEdgeStencilID,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >& macroFaceStencilID,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >& macroCellStencilID,
    const EdgeDoFForm&                                                                  form )
{
   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      /////////////////
      // Macro-edges //
      /////////////////

      for ( const auto& it : storage->getEdges() )
      {
         const auto& edge = *it.second;
         WALBERLA_ASSERT_GREATER( edge.getNumNeighborCells(), 0 );

         for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
         {
            const auto& cell = *( storage->getCell( edge.neighborCells().at( neighborCellID ) ) );

            const uint_t cellLocalEdgeID = cell.getLocalEdgeID( edge.getID() );
            const auto   basisInCell     = algorithms::getMissingIntegersAscending< 2, 4 >(
                {cell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
                 cell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

            const auto edgeAssemblyIndexInCell = indexing::basisConversion(
                indexing::Index( 0, 0, 0 ), basisInCell, {0, 1, 2, 3}, levelinfo::num_microedges_per_edge( level ) );

            auto& edgeToEdgeStencilMemory = edge.getData( macroEdgeStencilID )->getData( level );

            for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
               const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                   edgedof::EdgeDoFOrientation::X, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
               edgeToEdgeStencilMemory[neighborCellID][convertedCenterOrientation][leafOrientation] =
                   P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                       edgeAssemblyIndexInCell, convertedCenterOrientation, leafOrientation, cell, level, form );
            }
         }
      }

      /////////////////
      // Macro-faces //
      /////////////////

      if ( level >= 1 )
      {
         for ( const auto& it : storage->getFaces() )
         {
            const auto& face = *it.second;

            WALBERLA_ASSERT_GREATER( face.getNumNeighborCells(), 0 );

            for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
            {
               const auto& cell = *( storage->getCell( face.neighborCells().at( neighborCellID ) ) );

               const std::array< edgedof::EdgeDoFOrientation, 3 > faceOrientations = {
                   edgedof::EdgeDoFOrientation::X, edgedof::EdgeDoFOrientation::Y, edgedof::EdgeDoFOrientation::XY};

               const uint_t                  localFaceID          = cell.getLocalFaceID( face.getID() );
               const std::array< uint_t, 4 > localVertexIDsAtCell = {
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ),
                   6 - cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ) -
                       cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ) -
                       cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 )};

               const std::map< edgedof::EdgeDoFOrientation, indexing::Index > edgeAssemblyIndexInFace = {
                   {edgedof::EdgeDoFOrientation::X, indexing::Index( 0, 1, 0 )},
                   {edgedof::EdgeDoFOrientation::XY, indexing::Index( 0, 0, 0 )},
                   {edgedof::EdgeDoFOrientation::Y, indexing::Index( 1, 0, 0 )},
               };

               auto& edgeToEdgeStencilMemory = face.getData( macroFaceStencilID )->getData( level );
               for ( const auto& centerOrientation : faceOrientations )
               {
                  for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
                  {
                     const auto edgeAssemblyIndexInCell =
                         indexing::basisConversion( edgeAssemblyIndexInFace.at( centerOrientation ),
                                                    localVertexIDsAtCell,
                                                    {0, 1, 2, 3},
                                                    levelinfo::num_microedges_per_edge( level ) );
                     const auto convertedCenterOrientation =
                         edgedof::convertEdgeDoFOrientationFaceToCell( centerOrientation,
                                                                       localVertexIDsAtCell.at( 0 ),
                                                                       localVertexIDsAtCell.at( 1 ),
                                                                       localVertexIDsAtCell.at( 2 ) );
                     edgeToEdgeStencilMemory[neighborCellID][convertedCenterOrientation][leafOrientation] =
                         P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                             edgeAssemblyIndexInCell, convertedCenterOrientation, leafOrientation, cell, level, form );
                  }
               }
            }
         }
      }

      /////////////////
      // Macro-cells //
      /////////////////

      if ( level >= 1 )
      {
         for ( const auto& it : storage->getCells() )
         {
            const auto& cell = *it.second;

            auto& edgeToEdgeStencilMemory = cell.getData( macroCellStencilID )->getData( level );
            for ( const auto& centerOrientation : edgedof::allEdgeDoFOrientations )
            {
               if ( level == 1 && centerOrientation != edgedof::EdgeDoFOrientation::XYZ )
               {
                  continue;
               }
               for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
               {
                  const auto edgeToEdgeStencilMap = P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                      edgedof::macrocell::getInnerIndexByOrientation( centerOrientation ),
                      centerOrientation,
                      leafOrientation,
                      cell,
                      level,
                      form );
                  for ( const auto& stencilIt : edgeToEdgeStencilMap )
                  {
                     edgeToEdgeStencilMemory[centerOrientation][leafOrientation][stencilIt.first] = stencilIt.second;
                  }
               }
            }
         }
      }
   }
}

typedef EdgeDoFOperator< P2FenicsForm< hyteg::fenics::NoAssemble, fenics::NoAssemble > > GenericEdgeDoFOperator;

} // namespace hyteg
