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

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFApply.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {

template < class VertexDoFToEdgeDoFForm >
class VertexDoFToEdgeDoFOperator final : public Operator< P1Function< real_t >, EdgeDoFFunction< real_t > >
{
 public:
   VertexDoFToEdgeDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );
   VertexDoFToEdgeDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                               size_t                                     minLevel,
                               size_t                                     maxLevel,
                               const VertexDoFToEdgeDoFForm&              form );
   ~VertexDoFToEdgeDoFOperator() = default;

   void apply( const P1Function< real_t >&      src,
               const EdgeDoFFunction< real_t >& dst,
               size_t                           level,
               DoFType                          flag,
               UpdateType                       updateType = Replace ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< idx_t >&                  src,
                  const EdgeDoFFunction< idx_t >&             dst,
                  size_t                                      level,
                  DoFType                                     flag ) const;

   /// since the Vertex does not own any EdgeDoFs only edge, face and cell are needed
   const PrimitiveDataID< StencilMemory< real_t >, Edge >& getEdgeStencilID() const { return edgeStencilID_; }
   const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& getEdgeStencil3DID() const
   {
      return edgeStencil3DID_;
   }
   const PrimitiveDataID< StencilMemory< real_t >, Face >& getFaceStencilID() const { return faceStencilID_; }
   const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face >& getFaceStencil3DID() const
   {
      return faceStencil3DID_;
   }
   const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell >& getCellStencilID() const
   {
      return cellStencilID_;
   }

 private:
   void assembleStencils();

   PrimitiveDataID< StencilMemory< real_t >, Edge >                                      edgeStencilID_;
   PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge > edgeStencil3DID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >                                      faceStencilID_;
   PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face > faceStencil3DID_;
   PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell > cellStencilID_;

   VertexDoFToEdgeDoFForm form_;
};

template < typename VertexDoFToEdgeDoFForm >
void assembleVertexToEdgeStencils(
    const std::shared_ptr< PrimitiveStorage >&                                                   storage,
    const uint_t&                                                                                minLevel,
    const uint_t&                                                                                maxLevel,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& macroEdgeStencilID,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face >& macroFaceStencilID,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell >& macroCellStencilID,
    VertexDoFToEdgeDoFForm&                                                                      form )
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
            form.setGeometryMap( cell.getGeometryMap() );

            const uint_t cellLocalEdgeID = cell.getLocalEdgeID( edge.getID() );
            const auto   basisInCell     = algorithms::getMissingIntegersAscending< 2, 4 >(
                { cell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
                        cell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 ) } );

            const auto edgeAssemblyIndexInCell = indexing::basisConversion(
                indexing::Index( 0, 0, 0 ), basisInCell, { 0, 1, 2, 3 }, levelinfo::num_microedges_per_edge( level ) );

            auto& vertexToEdgeStencilMemory = edge.getData( macroEdgeStencilID )->getData( level );
            {
               const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                   edgedof::EdgeDoFOrientation::X, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
               vertexToEdgeStencilMemory[neighborCellID][convertedCenterOrientation] =
                   P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                       edgeAssemblyIndexInCell, convertedCenterOrientation, cell, level, form );
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
               form.setGeometryMap( cell.getGeometryMap() );

               const std::array< edgedof::EdgeDoFOrientation, 3 > faceOrientations = {
                   edgedof::EdgeDoFOrientation::X, edgedof::EdgeDoFOrientation::Y, edgedof::EdgeDoFOrientation::XY };

               const uint_t                  localFaceID          = cell.getLocalFaceID( face.getID() );
               const std::array< uint_t, 4 > localVertexIDsAtCell = {
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ),
                   6 - cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ) -
                       cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ) -
                       cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ) };

               const std::map< edgedof::EdgeDoFOrientation, indexing::Index > edgeAssemblyIndexInFace = {
                   { edgedof::EdgeDoFOrientation::X, indexing::Index( 0, 1, 0 ) },
                   { edgedof::EdgeDoFOrientation::XY, indexing::Index( 0, 0, 0 ) },
                   { edgedof::EdgeDoFOrientation::Y, indexing::Index( 1, 0, 0 ) },
               };

               auto& vertexToEdgeStencilMemory = face.getData( macroFaceStencilID )->getData( level );
               for ( const auto& centerOrientation : faceOrientations )
               {
                  const auto edgeAssemblyIndexInCell = indexing::basisConversion( edgeAssemblyIndexInFace.at( centerOrientation ),
                                                                                  localVertexIDsAtCell,
                                                                                  { 0, 1, 2, 3 },
                                                                                  levelinfo::num_microedges_per_edge( level ) );

                  const auto convertedCenterOrientation =
                      edgedof::convertEdgeDoFOrientationFaceToCell( centerOrientation,
                                                                    localVertexIDsAtCell.at( 0 ),
                                                                    localVertexIDsAtCell.at( 1 ),
                                                                    localVertexIDsAtCell.at( 2 ) );
                  vertexToEdgeStencilMemory[neighborCellID][convertedCenterOrientation] =
                      P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                          edgeAssemblyIndexInCell, convertedCenterOrientation, cell, level, form );
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
            form.setGeometryMap( cell.getGeometryMap() );

            auto& vertexToEdgeStencilMemory = cell.getData( macroCellStencilID )->getData( level );
            for ( const auto& centerOrientation : edgedof::allEdgeDoFOrientations )
            {
               if ( level == 1 && centerOrientation != edgedof::EdgeDoFOrientation::XYZ )
               {
                  continue;
               }

               const auto vertexToEdgeStencilMap = P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                   edgedof::macrocell::getInnerIndexByOrientation( centerOrientation ), centerOrientation, cell, level, form );
               for ( const auto& stencilIt : vertexToEdgeStencilMap )
               {
                  vertexToEdgeStencilMemory[centerOrientation][stencilIt.first] = stencilIt.second;
               }
            }
         }
      }
   }
}

namespace VertexDoFToEdgeDoF {

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroEdgeVertexDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive );

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroFaceVertexDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive );

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroCellVertexDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive );

} // namespace VertexDoFToEdgeDoF

typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< hyteg::fenics::NoAssemble, fenics::NoAssemble > >
                                                                                        GenericVertexDoFToEdgeDoFOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise > > VertexToEdgeDivTxOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise > > VertexToEdgeDivTyOperator;

typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >
    P1ToP2DivTxVertexToEdgeOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >
    P1ToP2DivTyVertexToEdgeOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >
    P1ToP2DivTzVertexToEdgeOperator;

} // namespace hyteg
